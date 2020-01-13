/**
 * @file    checker.hpp
 * @author  Vadim Piven <vadim@piven.tech>
 * @license Free use of this library is permitted under the
 * guidelines and in accordance with the MIT License (MIT).
 * @url     https://github.com/irreducible-polynoms/irrpoly
 */

#ifndef IRRPOLY_CHECKER_HPP
#define IRRPOLY_CHECKER_HPP

#include <thread>
#include <cassert>
#include <functional>

#ifdef PTHREAD

#include <pthread.h>

#else

#include <mutex>
#include <condition_variable>

#endif

namespace irrpoly {

    namespace detail {

        /// Класс для управления мьютексом и условной переменной.
        class sync {
        private:
#ifdef PTHREAD
            pthread_mutex_t mutex;
            pthread_cond_t cond;
#else
            ::std::mutex mutex;
            ::std::unique_lock<::std::mutex> lk;
            ::std::condition_variable cond;
#endif

        public:
            /// Инициализация и блокировка мьютекса.
            sync() noexcept(false) : cond(), mutex() {
#ifdef PTHREAD
                pthread_mutexattr_t attr;
                assert(!pthread_mutexattr_init(&attr));
                assert(!pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_ERRORCHECK));
                assert(!pthread_mutex_init(&mutex, &attr));
                assert(!pthread_cond_init(&cond, nullptr));
                assert(!pthread_mutexattr_destroy(&attr));
                pthread_mutex_lock(&mutex);
#else
                lk = ::std::unique_lock<::std::mutex>(mutex);
#endif
            }

            /// Освобождение ресурсов, здесь происходит разблокировка мьютекса.
            virtual
#ifndef PTHREAD
            ~sync() = default;

#else
            ~sync() noexcept(false) {
                pthread_mutex_unlock(&mutex);
                assert(!pthread_mutex_destroy(&mutex));
                assert(!pthread_cond_destroy(&cond));
            }

#endif

            /// Блокировка мьютекса.
            void lock() noexcept {
#ifdef PTHREAD
                pthread_mutex_lock(&mutex);
#else
                mutex.lock();
#endif
            }

            /// Разблокировка мьютекса.
            void unlock() noexcept {
#ifdef PTHREAD
                pthread_mutex_unlock(&mutex);
#else
                mutex.unlock();
#endif
            }

            /// Ожидание уведомления условной переменной.
            void wait() noexcept {
#ifdef PTHREAD
                pthread_cond_wait(&cond, &mutex);
#else
                cond.wait(lk);
#endif
            }

            /// Уведомление условной переменной.
            void signal() noexcept {
#ifdef PTHREAD
                pthread_cond_signal(&cond);
#else
                cond.notify_one();
#endif
            }
        };

    } // namespace detail

    /// Выполняет проверку на неприводимось и примитивность заданного многочлена над полем GF[P].
    template<typename value_type, typename result_type>
    class checker {
    public:
        /// Вид функции, генерирующей многочлены для проверки.
        typedef ::std::function<value_type()> input_func;

        /// Вид функции, выполняющей проверку и сохраняющей результат.
        typedef ::std::function<void(const value_type &, result_type &)> check_func;

        /// Вид функции, обрабатывающей результат проверки (если многочлен удовлетворяет
        /// требуемым условиям возвращает true, иначе - false).
        typedef ::std::function<bool(const value_type &, const result_type &)> callback_func;

    private:
        /// Потоки, непосредственно выполняющие проверку многочленов на неприводимость и примитивность.
        class node : public detail::sync {
            detail::sync *m; ///< Основной объект синхронизации
            check_func cf; ///< Основная функция проверки

            value_type val; ///< Входные данные
            result_type res; ///< Результат проверки

            bool _busy; ///< Поток занят полезной работой
            bool _terminate; ///< Поток должен быть завершён

#ifdef PTHREAD
            pthread_t thread;
#else
            ::std::thread thread;
#endif

            /**
             * Собственно функция, выполняющая проверку многочлена на неприводимость и примитивность.
             * Поток бесконечно ожидает получения новых данных. Если данные получены - выполняется их
             * проверка, выставление результата и уведомление условной переменной.
             * Кроме того, постоянно проверяется, не должен ли поток завершить работу.
             * Выход из функции прекращает работу потока и освобожает его ресурсы.
             * В случае с pthread для этого требуется выполнить pthread_exit.
             */
#ifdef PTHREAD

            static
            void *check(void *arg) noexcept {
                auto *sl = static_cast<node *>(arg);
#else
                void check() noexcept {
                    auto *sl = this;
#endif
                while (!sl->_terminate) {
                    while (!sl->_busy && !sl->_terminate) { sl->wait(); }
                    if (sl->_terminate) { break; }

                    sl->cf(sl->val, sl->res);

                    sl->m->lock();
                    sl->_busy = false;
                    sl->m->signal();
                    sl->m->unlock();
                }
#ifdef PTHREAD
                pthread_exit(nullptr);
#endif
            }

        public:
            explicit
            node(detail::sync *m) noexcept :
                    m(m), cf(), val(), res(), _busy(false), _terminate(false), thread() {
#ifdef PTHREAD
                pthread_attr_t thread_attr;
                assert(!pthread_attr_init(&thread_attr));
                assert(!pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_DETACHED));
                pthread_create(&thread, &thread_attr, &node::check, this);
                assert(!pthread_attr_destroy(&thread_attr));
#else
                thread = ::std::thread(&node::check, this);
                thread.detach();
#endif
            }

            /// Устанавливает функцию, которая будет вызываться в процессе проверки
            void set_check(check_func c) noexcept {
                cf = c;
            }

            /// Многочлен, проверка которого выполнялась.
            [[nodiscard]]
            const value_type &get() const {
                return val;
            }

            /**
             * Установка нового многочлена для проверки, сбрасывает результаты
             * предыдущей проверки и выставляет busy = true.
             */
            void set(value_type v) {
                lock();
                val = v;
                _busy = true;
                signal();
                unlock();
            }

            /// Возвращает текущее состояние потока.
            [[nodiscard]]
            bool busy() const noexcept {
                return _busy;
            }

            /// Устанавливает флаг, требующий завершить работу потока по завершении вычислений.
            void terminate() noexcept {
                lock();
                _terminate = true;
                signal();
                unlock();
            }

            /// Возвращает результат проверки текущего многочлена на неприводимость и примитивность.
            [[nodiscard]]
            const result_type &result() const noexcept {
                return res;
            }
        }; // class node

        detail::sync m; // master
        ::std::vector<node *> s; // slave

        /// Считает, сколько потоков заняты.
        unsigned countBusy() {
            unsigned res = 0;
            for (const auto *n : s) { res += n->busy(); }
            return res;
        }

    public:
        explicit
        checker(const unsigned n = ::std::thread::hardware_concurrency() - 1) noexcept : m(), s() {
            s.reserve(n);
            for (unsigned i = 0; i < n; ++i) {
                s.emplace_back(new node(&m));
            }
        }

        /// Основной цикл разделения работы на потоки.
        void check(uint64_t n, input_func in, check_func cf, callback_func back, const bool strict = true) noexcept {
            // заряжаем многочлены на проверку
            for (auto *sl : s) {
                sl->set_check(cf);
                sl->set(in());
            }
            while (n) {
                // ждём свободный поток
                while (countBusy() == s.size()) { m.wait(); }
                // находим свободные потоки и заряжаем новыми входными данными
                for (unsigned i = 0; i < s.size() && n; ++i) {
                    if (s[i]->busy()) { continue; }
                    if (back(s[i]->get(), s[i]->result())) { --n; }
                    s[i]->set(in());
                }
            }
            // ожидаем завершения всех потоков
            while (countBusy()) { m.wait(); }
            if (!strict) {
                // обрабатываем все проверенные многочлены, даже если их больше, чем требовалось найти
                for (auto *sl : s) { back(sl->get(), sl->result()); }
            }
        }

        /// Завершаем работу всех потоков.
        ~checker() noexcept {
            for (auto *sl : s) { sl->terminate(); }
        }
    };

} // namespace irrpoly

#endif //IRRPOLY_CHECKER_HPP
