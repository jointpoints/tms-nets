/**
 * @file    checker.hpp
 * @author  Vadim Piven <vadim@piven.tech>
 * @license Free use of this library is permitted under the
 * guidelines and in accordance with the MIT License (MIT).
 * @url     https://github.com/irreducible-polynoms/irrpoly
 */

#pragma once

#include <thread>
#include <cassert>
#include <functional>
#include <mutex>
#include <condition_variable>
#include <utility>
#include <vector>
#include <memory>
#include <optional>
#include <algorithm>

namespace irrpoly::multithread {

/**
 * Chains the input_fn, payload_fn and callback_fn provided. Executes this chain in several
 * threads at the same time.
 * TODO: change architecture to make checks faster.
 */
template<typename input_t, typename output_t>
class pipeline {
public:
    using input_fn = std::function<input_t()>;

    using payload_fn =
    std::function<void(const input_t &, std::optional<output_t> &)>;

    using callback_fn =
    std::function<bool(const input_t &, const output_t &)>;

private:
    class pod {
        std::shared_ptr<std::mutex> s_mutex;
        std::shared_ptr<std::condition_variable> s_cond;

        payload_fn m_pl;

        volatile bool m_terminate; ///< set to true when pod should be terminated
        std::optional<input_t> m_val;

        bool m_busy; ///< set to true during payload_fn execution, otherwise false
        std::optional<output_t> m_res;

        std::thread m_thread;
        std::mutex m_mutex;
        std::condition_variable m_cond;

        void execute() {
            std::unique_lock<std::mutex> lk(m_mutex);

            while (!m_terminate) {
                if (!(m_busy || m_terminate)) {
                    m_cond.wait(lk);
                    continue;
                }

                m_pl(m_val.value(), m_res);

                std::lock_guard<std::mutex> lg(*s_mutex);
                m_busy = false;
                s_cond->notify_one();
            }
        }

    public:
        pod(std::shared_ptr<std::mutex> s_mutex,
            std::shared_ptr<std::condition_variable> s_cond) :
            s_mutex(std::move(s_mutex)), s_cond(std::move(s_cond)),
            m_busy(false), m_terminate(false),
            m_val(), m_res() {
            m_thread = std::thread(&pod::execute, std::ref(*this));
        }

        ~pod() {
            // waits for thread to terminate, only then object resources will be freed
            m_thread.join();
        }

        void set_payload(payload_fn pl) {
            m_pl = pl;
        }

        [[nodiscard]]
        auto input() const -> const std::optional<input_t> & {
            return m_val;
        }

        void input(input_t v) {
            std::lock_guard<std::mutex> lg(m_mutex);
            m_val.emplace(std::move(v));
            m_res.reset();
            m_busy = true;
            m_cond.notify_one();
        }

        [[nodiscard]]
        auto is_busy() const -> bool {
            return m_busy;
        }

        void terminate() {
            std::lock_guard<std::mutex> lg(m_mutex);
            m_terminate = true;
            m_cond.notify_one();
        }

        [[nodiscard]]
        auto output() const -> const std::optional<output_t> & {
            return m_res;
        }

        void clean() {
            std::lock_guard<std::mutex> lg(m_mutex);
            m_val.reset();
            m_res.reset();
            m_busy = false;
        }
    }; // class node

    std::shared_ptr<std::mutex> s_mutex;
    std::shared_ptr<std::condition_variable> s_cond;

    std::vector<std::unique_ptr<pod>> m_pods;

    auto count_busy() -> unsigned {
        return std::count_if(m_pods.begin(), m_pods.end(),
                             std::mem_fn(&pod::is_busy));
    }

public:
    explicit
    pipeline(unsigned n = std::thread::hardware_concurrency()) {
        s_mutex = std::make_shared<std::mutex>();
        s_cond = std::make_shared<std::condition_variable>();

        if (n > 1) {
            m_pods.reserve(--n);
            for (unsigned i = 0; i < n; ++i) {
                m_pods.push_back(std::make_unique<pod>(s_mutex, s_cond));
            }
        }
    }

    void chain(input_fn in, payload_fn pl, callback_fn bk,
               const bool strict = true) {
        // if multithreading is unavailable - perform everything in main thread
        if (m_pods.empty()) {
            while (true) {
                auto input = in();
                std::optional<output_t> result;
                pl(input, result);
                if (bk(input, result.value())) {
                    return;
                }
            }
        }

        std::unique_lock<std::mutex> lk(*s_mutex);
        std::this_thread::yield();

        for (const auto &sl : m_pods) {
            sl->set_payload(pl);
            sl->input(in());
        }
        while (true) {
            while (count_busy() == m_pods.size()) {
                s_cond->wait(lk);
            }
            for (unsigned i = 0; i < m_pods.size(); ++i) {
                if (m_pods[i]->is_busy()) {
                    continue;
                }
                if (bk(m_pods[i]->input().value(),
                       m_pods[i]->output().value())) {
                    m_pods[i]->clean();
                    goto end;
                }
                m_pods[i]->input(in());
            }
        }
        end:
        while (count_busy()) {
            s_cond->wait(lk);
        }
        if (!strict) {
            // collect all results, by default excess results are discarded
            for (const auto &sl : m_pods) {
                if (sl->input().has_value() &&
                    sl->output().has_value()) {
                    bk(sl->input().value(), sl->output().value());
                    sl->clean();
                }
            }
        }
    }

    ~pipeline() {
        for (const auto &sl : m_pods) {
            sl->terminate();
        }
    }
};

} // namespace irrpoly::multithread
