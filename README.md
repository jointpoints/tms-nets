ИНСТРУКЦИЯ
==========




## Содержание

  * [Описание](#описание)
  * [Ветви репозитория](#ветви-репозитория)
  * [Технические требования](#технические-требования)
  * [Загрузка](#загрузка)
    * [Через git](#через-git)
    * [Через ZIP-архив](#через-zip-архив)
  * [Использование генератора](#использование-генератора)
    * [Краткое описание класса](#краткое-описание-класса)
    * [Включение файлов в Вашу программу](#включение-файлов-в-вашу-программу)
    * [Примеры](#примеры)
  * [Тестирование](#тестирование)
    * [Краткое описание тестов](#краткое-описание-тестов)
    * [Пределы применимости тестов](#пределы-применимости-тестов)
    * [Автоматическое тестирование](#автоматическое-тестирование)
    * [Сборка автоматического тестера с использованием CMake](#сборка-автоматического-тестера-с-использованием-cmake)
    * [Включение тестов в Вашу программу](#включение-тестов-в-вашу-программу)
    * [Уровни вербозности](#уровни-вербозности)
  * [Development Guides](#development-guides)
    * [Документация](#документация)
    * [Рекомендуемая литература](#рекомендуемая-литература)




## Описание

Данный репозиторий содержит в себе реализацию генератора участков (t, s)-последовательностей с основанием 2, созданного на основе
алгоритма, предложенного Гаральдом Нидеррайтером в 1987 году в работе "Low-Discrepancy and Low-Dispersion Sequences".




## Ветви репозитория

  * **master** - главная ветвь; в ней располагается последняя стабильная версия проекта;
  * **development** - ветвь, содержащая самые актуальные наработки;
  * **gh-pages** - техническая ветвь, хранящая в себе файлы документации.

Объединение **gh-pages** с другими ветвями, а также самостоятельное изменение её содержимого не допускается.




## Технические требования

Для компиляции и сборки программы, использующей данный генератор, необходимо иметь:

  * компилятор, поддерживающий C++17;
  * 64-битный процессор с 2 и более ядрами.




## Загрузка


#### Через git

Исполните команду

    git clone --single-branch --branch master https://github.com/jointpoints/tms-nets --recurse-submodules


#### Через ZIP-архив

  1. Скачайте архив 1 (кнопка **Clone or download**) [отсюда](https://github.com/jointpoints/tms-nets);
  2. Скачайте архив 2 (кнопка **Clone or download**) [отсюда](https://github.com/irreducible-polynoms/irrpoly);
  3. Распакуйте архив 1 в некоторую папку, например, `nets`;
  4. Распакуйте архив 2 в папку `nets/irrpoly`.




## Использование генератора


#### Краткое описание класса

Генератор последовательностей представлен шаблонным классом `sequences::Niederreiter<typename UIntType, unsigned int NBITS>`, где

  * `UIntType` - тип целочисленных переменных, который следует использовать при генерации;
  * `NBITS` - число доступных бит в `UIntType` для хранения промежуточных результатов расчётов, **максимальное значение** - `sizeof(UIntType) * 8`;

Сигнатуры конструкторов:

1.

    Niederreiter(uint32_t const dim, bool const in_parallel)

  * `dim` - размерность пространства;
  * `in_parallel` - указывает, следует ли генерировать неприводимые многочлены последовательно (`false`) или параллельно (`true`); **значение по умолчанию** - `false`.

Этот конструктор следует использовать, если требуется сгенерировать максимально плотно распределённые точки.

2.

    Niederreiter(std::vector<uint32_t> const &degrees_of_irred)

  * `degrees_of_irred` - вектор, задающий степени, которые следует использовать при генерации неприводимых многочленов; для *i*-й компоненты будет сгенерирован многочлен `degrees_of_irred[i-1]`-й степени при *i* от 1 до `degrees_of_irred.size()`.

Этот конструктор следует использовать, если требуется вручную настроить плотность по каждой компоненте.

3.

    Niederreiter(std::vector<irrpoly::polynomialgf<2>> const &polynomials)

  * `polynomials` - вектор неприводимых многочленов, которые следует использовать при генерации (t, s)-последовательности, где s = `polynomials.size()`; для *i*-й компоненты будет использован многочлен `polynomials[i-1]` при *i* от 1 до `polynomials.size()`.

Этот конструктор следует использовать, если требуется вручную настроить плотность по каждой компоненте наиболее свободно.


#### Включение файлов в Вашу программу

При написании программ, использующих данный генератор, потребуются файлы:

  * `niederreiter2.hpp`
  * `irrpoly/checker.hpp`
  * `irrpoly/gf.hpp`
  * `irrpoly/polynomial.hpp`
  * `irrpoly/polynomialgf.hpp`

В своём исходном коде пропишите строку

    #include "<your/path/to/our/files/>niederreiter2.hpp"


#### Примеры

В скором времени будут добавлены примеры использования с исчерпывающими комментариями.




## Тестирование


#### Краткое описание тестов

На данный момент реализованы тесты:

1.

    const bool niederreiter_check_definition(sequences::Niederreiter<UIntType, NBITS> *generator, uint8_t m)

Этот тест проверяет, являются ли 2^`m` точек, рассчитанных генератором `generator`, (t, m, s)-сетью в классическом смысле Нидеррайтера. Тест может давать ложноположительные результаты, если генерируемые точки не уникальны.

2.

    const bool niederreiter_check_uniqueness(sequences::Niederreiter<UIntType, NBITS> *generator, uint8_t m)

Этот тест проверяет, являются ли 2^`m` точек, рассчитанных генератором `generator`, покомпонентно уникальными.

Все тесты возвращают `true` при успешном прохождении и `false` при неуспешном.


#### Пределы применимости тестов

Реализованные тесты оптимально расходают ресурсы памяти компьютера, однако несмотря на это существуют некоторые ограничения:

  * В общем случае не рекомендуется при тестировании использовать значения `NBITS`, превышающие `32`;
  * Тест `niederreiter_check_definition` работает только при значениях `NBITS`, не превышающих `63`.


#### Автоматическое тестирование

В скором времени будет добавлена программа, порождающая случайные тестовые случаи и проводящая их контроль посредством всех реализованных тестов. На данный момент три тестовых случая находятся в `tests/main.cpp`.


#### Сборка автоматического тестера с использованием CMake

Для сборки автоматического тестера (или пока что замещающих его трёх демонстрационных тестовых случаев) можно воспользоваться утилитой CMake версии 3.8 или выше. Для сборки необходимо в командной строке перейти в папку с проектом и выполнить следующие команды:

    mkdir build
    cd build
    cmake .. -G "MinGW Makefiles"
    mingw32-make

После успешного выполнения в папке `build/test` будет находиться `test.exe`. Обратите внимание на то, что данный пример подразумевает использование компилятора MinGWw64. При использовании иного компилятора впишите после ключа `-G` вместо `MinGW Makefiles` формат Make-файлов, соответствующий Вашему компилятору.


#### Включение тестов в Вашу программу

Для собственного использования тестов понадобятся файлы:

  * `niederreiter2.hpp`
  * `irrpoly/checker.hpp`
  * `irrpoly/gf.hpp`
  * `irrpoly/polynomial.hpp`
  * `irrpoly/polynomialgf.hpp`
  * `tests/tests.hpp`
  * `tests/tests_routines.hpp`
  * `tests/tests_routines.cpp`

В своём исходном коде вместо

    #include "<your/path/to/our/files/>niederreiter2.hpp"

пропишите строку

    #include "<your/path/to/our/files/>tests/tests.hpp"


#### Уровни вербозности

Все реализованные тесты поддерживают пять уровней вербозности:

  * `0` - Вывод не производится.
  * `1` - Выводится имя теста и результат его выполнения.
  * `2` - Выводится имя теста, результат выполнения, а также расчётные параметры (t, m, s)-сети (для `niederreiter_check_definition`), либо число покомпонентно уникальных точек (для `niederreiter_check_uniqueness`).
  * `3` - Выводится детальная информация о результатах теста, включая все найденные ошибки в работе генератора.
  * `4` - Выводится детальная информация и о работе, и о результатах теста.

Задать желаемый уровень вербозности можно через макрос `VERBOSITY_LEVEL`, определить который нужно до включения файла `tests/tests.hpp`. По умолчанию уровень вербозности равен нулю.




## Development Guides


#### Документация

Документация к стабильной версии программы доступна [здесь](https://jointpoints.github.io/tms-nets/) (на английском языке),
генерируется автоматически с помощью утилиты Doxygen при каждом изменении в ветви **master**. Для полного понимания номенклатуры и алгоритма рекомендуется ознакомиться с литературой из списка ниже.


#### Рекомендуемая литература

  1. Harald Niederreiter "Low-Discrepancy and Low-Dispersion Sequences"
  2. Harald Niederreiter "Random Number Generation and Quasi-Monte Carlo Methods (особенно глава 4)
  3. Paul Brately, Bennett L. Fox, Harald Niederreiter "Implementation and Tests of Low-Discrepancy Sequences"
