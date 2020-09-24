MANUAL
======
[Э! А по-русски?!](#инструкция)




## Contents

  * [Introduction](#introduction)
  * [Repository branches](#repository-branches)
  * [Technical requirements](#technical-requirements)
  * [Download](#download)
    * [Via git](#via-git)
    * [As ZIP-archive](#as-zip-archive)
  * [The use of generator](#the-use-of-generator)
    * [Description of library structure](#description-of-library-structure)
    * [Description of `Niederreiter` class](#description-of-niederreiter-class)
    * [File inclusion into your program](#file-inclusion-into-your-program)
    * [Examples](#examples)
      * [Example 1. Basic usage](#example-1-basic-usage)
      * [Example 2. Replacement of the old generator](#example-2-replacement-of-the-old-generator)
      * [Example 3. Usage of handlers](#example-3-usage-of-handlers)
  * [Testing](#testing)
  * [Further development](#further-development)
    * [Documentation](#documentation)
    * [Recommended literature](#recommended-literature)




## Introduction

This repository contains a generator of digital (t, m, s)-nets in base 2 made on a basis of algorithm proposed by Harald Niederreiter in 1987 in his article "Low-Discrepancy and Low-Dispersion Sequences". Broadly speaking, (t, m, s)-nets are low-discrepancy discrete sets of points in *s*-dimensional unit cube. You may find more information about them in the sources listed in [Recommended literature](#recommended-literature).

[^ to the top ^](#contents)




## Repository branches

  * **master** — main branch; it contains the latest stable version of the library;
  * **development** — branch with the newest changes;
  * **gh-pages** — technical branch containing the documentation files;
  * **knowledge** — branch with the theoretical resources on the topic of (t, m, s)-nets made by the authors of this repository (in Russian).

Merge of **gh-pages** with any other branches or manual change of its contents is unacceptable.

[^ to the top ^](#contents)




## Technical requirements

This library binds you to have

  * C++17 compiler;
  * 64-bit processor with 2 or more cores.

[^ to the top ^](#contents)




## Download


#### Via git

Execute the following command line

    git clone --single-branch --branch master https://github.com/jointpoints/tms-nets

[^ to the top ^](#contents)


#### As ZIP-archive

  1. Download archive (**Code** button) from [here](https://github.com/jointpoints/tms-nets);
  2. Extract the archive into some folder.

[^ to the top ^](#contents)




## The use of generator


#### Description of library structure

This library is self-contained with its functionality being wholly provided within a sole namespace `tms`. Users are strongly advised to be informed about the following components of namespace `tms`:

  * `Niederreiter` — class of generator itself, read the section below to learn more;
  * `Real` — type of floating-point number;
  * `Point` — type of *s*-dimensional point with components of type `Real`;
  * `Polynomial` — class of polynomials over **F₂**, read its [source](https://github.com/irreducible-polynoms/irrpoly/) to learn more.

[^ to the top ^](#contents)


#### Description of `Niederreiter` class

The generator of sequences is represented by a template class `tms::Niederreiter<typename UIntType>`, where

  * `UIntType` — unsigned integral type, that is used during the process of generation for the storage of temporary data.

Constructors signatures:

1.

    Niederreiter(BasicInt const nbits, BasicInt const dim, bool const in_parallel)

  * `nbits` — *m* parameter of the net (binary logarithm of number of points in a desired net), **maximum value** — `sizeof(UIntType) * 8`;
  * `dim` — *s* parameter of the net (dimension of the unit cube to be filled with points);
  * `in_parallel` — specifies whether the irreducible polynomials should be generated consecutively (`false`) or concurrently (`true`); **default value** — `false`.

Constructs the generator of (t, m, s)-net with specified values of *m*, *s*, and with induced least possible value of *t*.

2.

    Niederreiter(BasicInt const nbits, std::vector<BasicInt> const &degrees_of_irrpolys, std::vector< std::vector<uintmax_t> > const &matrix_of_initial_values)
    Niederreiter(BasicInt const nbits, std::initializer_list<BasicInt> const &degrees_of_irrpolys, std::initializer_list< std::vector<uintmax_t> > const &matrix_of_initial_values)

  * `nbits` — *m* parameter of the net, **maximum value** — `sizeof(UIntType) * 8`;
  * `degrees_of_irrpolys` — vector of degrees of irreducible polynomials that should be used for generation of points; for *i*-th component of `degrees_of_irrpolys.size()`-dimensional space a polynomial of degree `degrees_of_irrpolys[i-1]` will be generated automatically;
  * `matrix_of_initial_values` — matrix of initial values for all recursive sequences, **default value** — empty matrix.

Constructs the generator of (t, m, s)-net with specified *m* parameter, values of *s* degrees of initial irreducible polynomials,  and (optional, for advanced users) initial values of all recursive sequences, defining the generation matrices.

3.

    Niederreiter(BasicInt const nbits, std::vector< std::vector<uintmax_t> > const &irrpolys_coeffs, std::vector< std::vector<uintmax_t> > const &matrix_of_initial_values)
    Niederreiter(BasicInt const nbits, std::initializer_list< std::vector<uintmax_t> > const &irrpolys_coeffs, std::initializer_list< std::vector<uintmax_t> > const &matrix_of_initial_values)

  * `nbits` — *m* parameter of the net, **maximum value** — `sizeof(UIntType) * 8`;
  * `irrpolys_coeffs` — coefficients of irreducible polynomials over **F₂** which should be used during the construction of (t, m, s)-net; for *i*-th component of `irrpolys_coeffs.size()`-dimensional space the `irrpolys_coeffs[i-1]`-th polynomial will be used;
  * `matrix_of_initial_values` — matrix of initial values for all recursive sequences, **default value** — empty matrix.

Constructs the generator of (t, m, s)-net with specified *m* parameter, *s* initial irreducible polynomials, and (optional, for advanced users) initial values of all recursive sequences, defining the generation matrices.

[^ to the top ^](#contents)


#### File inclusion into your program

To use our generator include the `include/tms-nets/niederreiter2.hpp` file.

[^ to the top ^](#contents)


#### Examples

###### Example 1. Basic usage

The code snippet below generates and prints all 16 points of (t, 4, 3)-net where optimal *t* is selected automatically.

    uint32_t s = 3;
    uint32_t m = 4; // m = log₂ (amount of points in net) = log₂ 16
    tms::Niederreiter<uint64_t> generator(m, s); // constructor #1
    tms::Point                  tmp_point;
    for (uint32_t point_i = 0; point_i < (1U << m); ++point_i)
    {
        tmp_point = generator.generate_point(point_i); // generates point_i-th point of the net
        for (uint32_t dim_i = 0; dim_i < s; ++dim_i)
            std::cout << tmp_point[dim_i] << '\t';
        std::cout << '\n';
    }

Output:

	0       0       0
	0.5     0.5     0.25
	0.75    0.25    0.5
	0.25    0.75    0.75
	0.375   0.375   0.3125
	0.875   0.875   0.0625
	0.625   0.125   0.8125
	0.125   0.625   0.5625
	0.1875  0.3125  0.9375
	0.6875  0.8125  0.6875
	0.9375  0.0625  0.4375
	0.4375  0.5625  0.1875
	0.3125  0.1875  0.625
	0.8125  0.6875  0.875
	0.5625  0.4375  0.125
	0.0625  0.9375  0.375

###### Example 2. Replacement of the old generator

The library provides solution for cases when developers need to use multiple generators with different parametes consecutively. This might be achieved with the assignment operator `=` — generator that was kept in the assigned variable before is destroyed. E.g., the code snippet below generates 5 points from 5 different (t, 32, 4)-nets.

    uint32_t                    amount   = 5;
    std::vector<tms::BasicInt>  poly_deg = {1, 2, 3, 4};
    tms::Niederreiter<uint32_t> generator(32, poly_deg); // constructor #2
    tms::Point                  tmp_point;
    for (uint32_t point_i = 0; point_i < amount; ++point_i)
    {
        tmp_point = generator.generate_point(point_i);
        for (uint32_t dim_i = 0; dim_i < poly_deg.size(); ++dim_i)
            std::cout << tmp_point[dim_i] << '\t';
        std::cout << '\n';
        poly_deg = {point_i+2, point_i+3, point_i+4, point_i+5}; // provide new degrees of polynomials
        generator = tms::Niederreiter<uint32_t>(32, poly_deg); // old generator is destroyed, replaced by new one
    }

Output:

    0       0       0        0
    0.25    0.125   0.0625   0.03125
    0.375   0.1875  0.09375  0.046875
    0.125   0.0625  0.03125  0.015625
    0.1875  0.09375 0.046875 0.0234375

###### Example 3. Usage of handlers

The library allows to set user-defined functions to handle each generated point of sequence. For instance, the code snippet below prints the first coordinare of each point.

    // Handlers must have the following signature
    void my_point_handler(tms::Point const &new_point, tms::CountInt point_i)
    {
        std::cout << "point[" << point_i << "]'s first component is " << new_point[0] << '\n';
        return;
    }
    
    <...>
    
    uint32_t                    start_i  = 3;
    uint32_t                    amount   = 5;
    std::vector<tms::BasicInt>  poly_deg = {3, 2, 8, 5};
    tms::Niederreiter<uint32_t> generator(20, poly_deg); // constructor #2
    tms::Point                  tmp_point;
    generator.for_each_point(my_point_handler, amount, start_i); // generates (amount) points starting with (start_i)-th and sends them to (my_point_handler)

Output:

    point[3]'s first component is 0.25
    point[4]'s first component is 0.875
    point[5]'s first component is 0.75
    point[6]'s first component is 0.5
    point[7]'s first component is 0.625

[^ to the top ^](#contents)




## Testing

For information about testing see [TsTests Usage and Development Guide](https://github.com/jointpoints/tms-nets/blob/master/tests/README.md).

[^ to the top ^](#contents)




## Further development


#### Documentation

Documentation for a stable version may be found [here](https://jointpoints.github.io/tms-nets/).
It is generated automatically by Doxygen after every change in **master** branch. For the complete understanding of notation and algorithms it is highly recommended to get acquainted with the literature listed below.

[^ to the top ^](#contents)


#### Recommended literature

  1. Documents from **knowledge** branch (in Russian)
  1. Harald Niederreiter "Low-Discrepancy and Low-Dispersion Sequences"
  2. Harald Niederreiter "Random Number Generation and Quasi-Monte Carlo Methods" (especially chapter 4)
  3. Paul Brately, Bennett L. Fox, Harald Niederreiter "Implementation and Tests of Low-Discrepancy Sequences"

[^ to the top ^](#contents)







ИНСТРУКЦИЯ
==========
[English, please.](#manual)




## Содержание

  * [Введение](#введение)
  * [Ветви репозитория](#ветви-репозитория)
  * [Технические требования](#технические-требования)
  * [Загрузка](#загрузка)
    * [Через git](#через-git)
    * [Через ZIP-архив](#через-zip-архив)
  * [Использование генератора](#использование-генератора)
    * [Описание структуры библиотеки](#описание-структуры-библиотеки)
    * [Описание класса `Niederreiter`](#описание-класса-niederreiter)
    * [Включение файлов в Вашу программу](#включение-файлов-в-вашу-программу)
    * [Примеры](#примеры)
      * [Пример 1. Простейший способ использования](#пример-1-простейший-способ-использования)
      * [Пример 2. Замена старого генератора новым](#пример-2-замена-старого-генератора-новым)
      * [Пример 3. Использование обработчиков генерации](#пример-3-использование-обработчиков-генерации)
  * [Тестирование](#тестирование)
  * [Для дальнейшей разработки](#для-дальнейшей-разработки)
    * [Документация](#документация)
    * [Рекомендуемая литература](#рекомендуемая-литература)




## Введение

Данный репозиторий содержит в себе реализацию генератора цифровых (t, m, s)-сетей с основанием 2, созданного на основе алгоритма, предложенного Гаральдом Нидеррайтером в 1987 году в работе "Low-Discrepancy and Low-Dispersion Sequences". (t, m, s)-сети, вообще говоря, — это дискретные множества точек, однородно распределённых по *s*-мерному единичному кубу. Подробнее о них можно узнать в ресурсах, указанных в [рекомендуемой литературе](#рекомендуемая-литература).

[^ наверх ^](#содержание)




## Ветви репозитория

  * **master** — главная ветвь; в ней располагается последняя стабильная версия проекта;
  * **development** — ветвь, содержащая самые актуальные наработки;
  * **gh-pages** — техническая ветвь, хранящая в себе файлы документации;
  * **knowledge** — ветвь, содержащая теоретические источники по (t, m, s)-сетям, созданные авторами данного репозитория.

Объединение **gh-pages** с другими ветвями, а также самостоятельное изменение её содержимого не допускается.

[^ наверх ^](#содержание)




## Технические требования

Для работы библиотеки необходимо иметь:

  * компилятор, поддерживающий C++17;
  * 64-битный процессор с 2 и более ядрами.

[^ наверх ^](#содержание)




## Загрузка


#### Через git

Исполните команду

    git clone --single-branch --branch master https://github.com/jointpoints/tms-nets

[^ наверх ^](#содержание)


#### Через ZIP-архив

  1. Скачайте архив (кнопка **Code**) [отсюда](https://github.com/jointpoints/tms-nets);
  2. Распакуйте архив в некоторую папку.

[^ наверх ^](#содержание)




## Использование генератора


#### Описание структуры библиотеки

Функционал данной библиотеки заключён в единственное пространство имён `tms`, содержащее в себе все необходимые для работы определения. Пользователям рекомендуется знать о следующих компонентах пространства имён `tms`:

  * `Niederreiter` — непосредственно класс генератора, подробнее о его использовании написано ниже;
  * `Real` — тип вещественных чисел с плавающей точкой;
  * `Point` — тип точки в *s*-мерном пространстве с компонентами типа `Real`;
  * `Polynomial` — класс многочленов над полем **F₂**, подробнее о его использовании написано в [источнике](https://github.com/irreducible-polynoms/irrpoly/).

[^ наверх ^](#содержание)


#### Описание класса `Niederreiter`

Генератор последовательностей представлен шаблонным классом `tms::Niederreiter<typename UIntType>`, где

  * `UIntType` — тип целочисленных переменных, который следует использовать в процессе генерации для хранения промежуточных результатов расчётов.

Сигнатуры конструкторов:

1.

    Niederreiter(BasicInt const nbits, BasicInt const dim, bool const in_parallel)

  * `nbits` — параметр *m* сети (двоичный логарифм числа точек в желаемой сети), **максимальное значение** — `sizeof(UIntType) * 8`;
  * `dim` — параметр *s* сети (размерность единичного куба, заполняемого точками);
  * `in_parallel` — указывает, следует ли генерировать неприводимые многочлены последовательно (`false`) или параллельно (`true`); **значение по умолчанию** — `false`.

Порождает генератор (t, m, s)-сетей с указанными значениями *m*, *s* и с минимальным возможным *t*.

2.

    Niederreiter(BasicInt const nbits, std::vector<BasicInt> const &degrees_of_irrpolys, std::vector< std::vector<uintmax_t> > const &matrix_of_initial_values)
    Niederreiter(BasicInt const nbits, std::initializer_list<BasicInt> const &degrees_of_irrpolys, std::initializer_list< std::vector<uintmax_t> > const &matrix_of_initial_values)

  * `nbits` — параметр *m* сети, **максимальное значение** — `sizeof(UIntType) * 8`;
  * `degrees_of_irrpolys` — степени, которые следует использовать при генерации неприводимых многочленов; для *i*-й компоненты `degrees_of_irrpolys.size()`-мерного пространства будет сгенерирован многочлен `degrees_of_irrpolys[i-1]`-й степени.
  * `matrix_of_initial_values` — матрица инициализирующих значений для всех рекуррентных последовательностей, **значение по умолчанию** — пустая матрица.

Порождает генератор (t, m, s)-сетей с указанными значениями *m*, *s* степеней неприводимых многочленов и (опционально, для продвинутых пользователей) инициализирующих элементов всех рекуррентных последовательностей.

3.

    Niederreiter(BasicInt const nbits, std::vector< std::vector<uintmax_t> > const &irrpolys_coeffs, std::vector< std::vector<uintmax_t> > const &matrix_of_initial_values)
    Niederreiter(BasicInt const nbits, std::initializer_list< std::vector<uintmax_t> > const &irrpolys_coeffs, std::initializer_list< std::vector<uintmax_t> > const &matrix_of_initial_values)

  * `nbits` — параметр *m* сети, **максимальное значение** — `sizeof(UIntType) * 8`;
  * `irrpolys_coeffs` — коэффициенты неприводимых многочленов над **F₂**, которые следует использовать при генерации (t, m, s)-сети; для *i*-й компоненты `irrpolys_coeffs.size()`-мерного пространства будет использован `irrpolys_coeffs[i-1]`-ый многочлен;
  * `matrix_of_initial_values` — матрица инициализирующих значений для всех рекуррентных последовательностей, **значение по умолчанию** — пустая матрица.

Порождает генератор (t, m, s)-сетей с указанными значениями *m*, коэффициентов *s* неприводимых многочленов и (опционально, для продвинутых пользователей) инициализирующих элементов всех рекуррентных последовательностей.

[^ наверх ^](#содержание)


#### Включение файлов в Вашу программу

Для использования данного генератора подключите файл `include/tms-nets/niederreiter2.hpp`.

[^ наверх ^](#содержание)


#### Примеры

###### Пример 1. Простейший способ использования

Участок кода ниже генерирует и выводит все 10 точек (t, 4, 3)-сети, где оптимальное *t* выбирается автоматически.

    uint32_t s = 3;
    uint32_t m = 4; // m = log₂ (число точек в сети) = log₂ 16
    tms::Niederreiter<uint64_t> generator(m, s); // конструктор #1
    tms::Point                  tmp_point;
    for (uint32_t point_i = 0; point_i < (1U << m); ++point_i)
    {
        tmp_point = generator.generate_point(point_i); // генерирует point_i-ую точку сети
        for (uint32_t dim_i = 0; dim_i < s; ++dim_i)
            std::cout << tmp_point[dim_i] << '\t';
        std::cout << '\n';
    }

Вывод программы:

	0       0       0
	0.5     0.5     0.25
	0.75    0.25    0.5
	0.25    0.75    0.75
	0.375   0.375   0.3125
	0.875   0.875   0.0625
	0.625   0.125   0.8125
	0.125   0.625   0.5625
	0.1875  0.3125  0.9375
	0.6875  0.8125  0.6875
	0.9375  0.0625  0.4375
	0.4375  0.5625  0.1875
	0.3125  0.1875  0.625
	0.8125  0.6875  0.875
	0.5625  0.4375  0.125
	0.0625  0.9375  0.375

###### Пример 2. Замена старого генератора новым

Библиотека предусматривает случаи, когда разработчикам требуется использовать несколько генераторов с разными параметрами последовательно. Для этого можно воспользоваться оператором присваивания `=` — генератор, хранившийся в памяти до этого, уничтожается. Так, участок кода ниже генерирует 5 точек из 5 различных (t, 32, 4)-сетей.

    uint32_t                    amount   = 5;
    std::vector<tms::BasicInt>  poly_deg = {1, 2, 3, 4};
    tms::Niederreiter<uint32_t> generator(32, poly_deg); // конструктор #2
    tms::Point                  tmp_point;
    for (uint32_t point_i = 0; point_i < amount; ++point_i)
    {
        tmp_point = generator.generate_point(point_i);
        for (uint32_t dim_i = 0; dim_i < poly_deg.size(); ++dim_i)
            std::cout << tmp_point[dim_i] << '\t';
        std::cout << '\n';
        poly_deg = {point_i+2, point_i+3, point_i+4, point_i+5}; // новые степени многочленов
        generator = tms::Niederreiter<uint32_t>(32, poly_deg); // старый генератор уничтожается, замещается новым
    }

Вывод программы:

    0       0       0        0
    0.25    0.125   0.0625   0.03125
    0.375   0.1875  0.09375  0.046875
    0.125   0.0625  0.03125  0.015625
    0.1875  0.09375 0.046875 0.0234375

###### Пример 3. Использование обработчиков генерации

Библиотека позволяет назначать пользовательские функции, которые будут обрабатывать каждую сгенерированную точку из последовательности. Например, участок программы ниже выводит первую координату каждой точки.

    // Обработчики должны иметь следующую сигнатуру
    void my_point_handler(tms::Point const &new_point, tms::CountInt point_i)
    {
        std::cout << "point[" << point_i << "]'s first component is " << new_point[0] << '\n';
        return;
    }
    
    <...>
    
    uint32_t                    start_i  = 3;
    uint32_t                    amount   = 5;
    std::vector<tms::BasicInt>  poly_deg = {3, 2, 8, 5};
    tms::Niederreiter<uint32_t> generator(20, poly_deg); // конструктор #2
    tms::Point                  tmp_point;
    generator.for_each_point(my_point_handler, amount, start_i); // генерирует (amount) точек, начиная со (start_i)-й, и отправляет их в (my_point_handler)

Вывод программы:

    point[3]'s first component is 0.25
    point[4]'s first component is 0.875
    point[5]'s first component is 0.75
    point[6]'s first component is 0.5
    point[7]'s first component is 0.625

[^ наверх ^](#содержание)




## Тестирование

Информация о тестировании может быть найдена в документе [TsTests Usage and Development Guide](https://github.com/jointpoints/tms-nets/blob/master/tests/README.md) (на английском языке).

[^ наверх ^](#содержание)




## Для дальнейшей разработки


#### Документация

Документация к стабильной версии программы доступна [здесь](https://jointpoints.github.io/tms-nets/) (на английском языке),
генерируется автоматически с помощью утилиты Doxygen при каждом изменении в ветви **master**. Для полного понимания номенклатуры и алгоритма рекомендуется ознакомиться с литературой из списка ниже.

[^ наверх ^](#содержание)


#### Рекомендуемая литература

  1. Документы из ветви **knowledge**
  2. Harald Niederreiter "Low-Discrepancy and Low-Dispersion Sequences"
  3. Harald Niederreiter "Random Number Generation and Quasi-Monte Carlo Methods" (особенно глава 4)
  4. Paul Brately, Bennett L. Fox, Harald Niederreiter "Implementation and Tests of Low-Discrepancy Sequences"

[^ наверх ^](#содержание)
