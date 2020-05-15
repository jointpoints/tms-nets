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

This repository contains a generator of digital (t, s)-sequences in base 2 made on a basis of algorithm proposed by Harald Niederreiter in 1987 in his article "Low-Discrepancy and Low-Dispersion Sequences". (t, s)-sequences are a handy tool for the construction of (t, m, s)-nets — low-discrepancy discrete sets of points in *s*-dimensional unit cube. You may find more information about (t, s)-sequences and (t, m, s)-nets in the sources listed in [Recommended literature](#recommended-literature).

[^ to the top ^](#contents)




## Repository branches

  * **master** — main branch; it contains the latest stable version of the library;
  * **development** — branch with the newest changes;
  * **gh-pages** — technical branch containing the documentation files.

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

  1. Download archive (**Clone or download** button) from [here](https://github.com/jointpoints/tms-nets);
  2. Extract the archive into some folder.

[^ to the top ^](#contents)




## The use of generator


#### Description of library structure

The functionality of this library is wholly provided within a sole namespace `sequences` containing all necessary definitions. Users are strongly advised to be informed about the following components of namespace `sequences`:

  * `Niederreiter` — class of generator itself, read the section below to learn more;
  * `Real` — type of floating-point number;
  * `Point` — type of *s*-dimensional point with components of type `Real`;
  * `Polynom` — class of polynomials over **F₂**, read its [source](https://github.com/irreducible-polynoms/irrpoly/) to learn more.

[^ to the top ^](#contents)


#### Description of `Niederreiter` class

The generator of sequences is represented by a template class `sequences::Niederreiter<typename UIntType, unsigned int NBITS>`, where

  * `UIntType` — unsigned integral type, that is used during the process of generation for the storage of temporary data;
  * `NBITS` — number of bits in `UIntType` that can be used for the storage of temporary data, **maximum value** — `sizeof(UIntType) * 8`.

Digital (t, s)-sequences have period. In our implementation the period of sequence is `2^NBITS`, which is why the more unique points are needed to be generated, the greater value of `NBITS` and, hence, the more capacious `UIntType` are needed to be used.

Constructors signatures:

1.

    Niederreiter(BasicInt const dim, bool const in_parallel)

  * `dim` — dimension of the unit cube to be filled with points;
  * `in_parallel` — specifies whether the irreducible polynomials should be generated consecutively (`false`) or concurrently (`true`); **default value** — `false`.

This constructor allows to generate points with the lowest possible discrepancy.

2.

    Niederreiter(std::vector<BasicInt> const &degrees_of_irred)

  * `degrees_of_irred` — vector of degrees of irreducible polynomials that should be used for generation of points; for *i*-th component of `degrees_of_irred.size()`-dimensional space a polynomial of degree `degrees_of_irred[i-1]` will be generated automatically.

This constructor is to be used when the manual configuration of discrepancy along any of the axes is needed to be done. The greater the degree of the *i*-th polynomial, the greater the discrepancy over the *i*-th component will be. The first constructor minimises these degrees which helps to achieve the best results.

3.

    Niederreiter(std::vector<Polynom> const &polynomials)

  * `polynomials` — vector of irreducible polynomials over **F₂** which should be used during the construction of (t, s)-sequence; for *i*-th component of `polynomials.size()`-dimensional space the `polynomials[i-1]`-th polynomial will be used.

This constructor is to be used when the manual configuration of discrepancy along any of the axes is needed to be done. Unlike the second constructor, this constructor expects users to specify the irreducible polynomials explicitly.

[^ to the top ^](#contents)


#### File inclusion into your program

To write programs that use our generator the files from `include` folder will be needed.

Write in your code

    #include "<your/path/to/our/files/>include/niederreiter2.hpp"

[^ to the top ^](#contents)


#### Examples

###### Example 1. Basic usage

The code snippet below generates and prints the first 10 points of 3-dimensional (t, s)-sequence.

    uint32_t dim    = 3;
    uint32_t amount = 10;
    sequences::Niederreiter<uint64_t, 64>   generator(dim);
    sequences::Point                        tmp_point;
    for (uint32_t i = 0; i < amount; ++i)
    {
        tmp_point = generator.get_point_real(i); // generates i-th point of sequence
        for (uint32_t j = 0; j < dim; ++j)
            std::cout << tmp_point[j] << '\t';
        std::cout << '\n';
    }

Output:

    0       0       0
    0.5     0.5     0.75
    0.75    0.25    0.3125
    0.25    0.75    0.5625
    0.375   0.375   0.875
    0.875   0.875   0.125
    0.625   0.125   0.6875
    0.125   0.625   0.4375
    0.1875  0.3125  0.515625
    0.6875  0.8125  0.265625

###### Example 2. Replacement of the old generator

The library provides for cases when developers need to use multiple generators with different parametes consecutively. This might be achieved with the assignment operator `=` — generator that was kept in the assigned variable before is destroyed. E.g., the code snippet below generates 5 points from 5 different sequences in a 4-dimensional space.

    uint32_t                          amount   = 5;
    std::vector<sequences::BasicInt>  poly_deg = {1, 2, 3, 4};
    sequences::Niederreiter<uint32_t, 32>  generator(poly_deg);
    sequences::Point                       tmp_point;
    for (uint32_t i = 0; i < amount; ++i)
    {
        tmp_point = generator.get_point_real(i);
        for (uint32_t j = 0; j < poly_deg.size(); ++j)
            std::cout << tmp_point[j] << '\t';
        std::cout << '\n';
        poly_deg = {i+2, i+3, i+4, i+5}; // provide new degrees of polynomials
        generator = sequences::Niederreiter<uint32_t, 32>(poly_deg); // old generator is destroyed, replaced by new one
    }

Output:

    0               0               0               0
    0.75            0.875           0.9375          0.96875
    0.140625        0.0664062       0.0322266       0.0158691
    0.878906        0.938477        0.968994        0.984436
    0.0644531       0.0317383       0.0157471       0.00784302

###### Example 3. Usage of handlers

The library allows to set user-defined functions to handle each generated point of sequence. For instance, the code snippet below prints the first coordinare of each point.

    // Handlers must have the following signature
    void my_point_handler(sequences::Point const &new_point, sequences::CountInt point_i)
    {
        std::cout << "point[" << point_i << "]'s first component is " << new_point[0] << '\n';
        return;
    }
    
    <...>
    
    uint32_t                          start_i  = 3;
    uint32_t                          amount   = 5;
    std::vector<sequences::BasicInt>  poly_deg = {3, 2, 8, 5};
    sequences::Niederreiter<uint32_t, 32>  generator(poly_deg);
    sequences::Point                       tmp_point;
    generator.for_each_point_real(my_point_handler, amount, start_i); // Generates (amount) points starting with (start_i)-th and sends them to (my_point_handler)

Output:

    point[3]'s first component is 0.765625
    point[4]'s first component is 0.28125
    point[5]'s first component is 0.65625
    point[6]'s first component is 0.421875
    point[7]'s first component is 0.546875

[^ to the top ^](#contents)




## Testing

For information about testing see [TsTests Usage and Development Guide](https://github.com/jointpoints/tms-nets/blob/master/tests/README.md).

[^ to the top ^](#contents)




## Further development


#### Documentation

Documentation for a stable version may be found [here](https://jointpoints.github.io/tms-nets/) (partly in Russian).
It is generated automatically by Doxygen after every change in **master** branch. For the complete understanding of notation and algorithms it is highly recommended to get acquainted with the literature listed below.

[^ to the top ^](#contents)


#### Recommended literature

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

Данный репозиторий содержит в себе реализацию генератора цифровых (t, s)-последовательностей с основанием 2, созданного на основе алгоритма, предложенного Гаральдом Нидеррайтером в 1987 году в работе "Low-Discrepancy and Low-Dispersion Sequences". (t, s)-последовательности представляют собой удобный инструмент построения (t, m, s)-сетей — дискретного множества точек, однородно распределённых по *s*-мерному единичному кубу. Подробнее о (t, s)-последовательностях и (t, m, s)-сетях — в ресурсах, указанных в [рекомендуемой литературе](#рекомендуемая-литература).

[^ наверх ^](#содержание)




## Ветви репозитория

  * **master** — главная ветвь; в ней располагается последняя стабильная версия проекта;
  * **development** — ветвь, содержащая самые актуальные наработки;
  * **gh-pages** — техническая ветвь, хранящая в себе файлы документации.

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

  1. Скачайте архив (кнопка **Clone or download**) [отсюда](https://github.com/jointpoints/tms-nets);
  2. Распакуйте архив в некоторую папку.

[^ наверх ^](#содержание)




## Использование генератора


#### Описание структуры библиотеки

Функционал данной библиотеки заключён в единственное пространство имён `sequences`, содержащее в себе все необходимые для работы определения. Пользователям рекомендуется знать о следующих компонентах пространства имён `sequences`:

  * `Niederreiter` — непосредственно класс генератора, подробнее о его использовании написано ниже в этом разделе;
  * `Real` — тип вещественных чисел с плавающей точкой;
  * `Point` — тип точки в *s*-мерном пространстве с компонентами типа `Real`;
  * `Polynom` — класс многочленов над полем **F₂**, подробнее о его использовании написано в [источнике](https://github.com/irreducible-polynoms/irrpoly/).

[^ наверх ^](#содержание)


#### Описание класса `Niederreiter`

Генератор последовательностей представлен шаблонным классом `sequences::Niederreiter<typename UIntType, unsigned int NBITS>`, где

  * `UIntType` — тип целочисленных переменных, который следует использовать в процессе генерации для хранения промежуточных результатов расчётов;
  * `NBITS` — число доступных бит в `UIntType` для хранения промежуточных результатов расчётов, **максимальное значение** — `sizeof(UIntType) * 8`.

Цифровые (t, s)-последовательности периодичны. В нашей реализации период последовательности равен `2^NBITS`, поэтому чем больше уникальных точек необходимо сгенерировать, тем большее значение `NBITS` и, как следствие, тем более объёмный тип `UIntType` необходимо указывать.

Сигнатуры конструкторов:

1.

    Niederreiter(BasicInt const dim, bool const in_parallel)

  * `dim` — размерность пространства, единичный куб внутри которого необходимо заполнить;
  * `in_parallel` — указывает, следует ли генерировать неприводимые многочлены последовательно (`false`) или параллельно (`true`); **значение по умолчанию** — `false`.

Этот конструктор позволяет генерировать максимально однородно заполняющие пространство точки.

2.

    Niederreiter(std::vector<BasicInt> const &degrees_of_irred)

  * `degrees_of_irred` — вектор, задающий степени, которые следует использовать при генерации неприводимых многочленов; для *i*-й компоненты `degrees_of_irred.size()`-мерного пространства будет сгенерирован многочлен `degrees_of_irred[i-1]`-й степени.

Этот конструктор следует использовать, если требуется вручную настроить однородность точек вдоль каждой оси. Чем выше степень *i*-го многочлена, тем менее однородно по *i*-ой компоненте распределены точки. Первый конструктор минимизирует все используемые степени, что позволяет добиться наилучшего результата по данному показателю.

3.

    Niederreiter(std::vector<Polynom> const &polynomials)

  * `polynomials` — вектор неприводимых многочленов над **F₂**, которые следует использовать при генерации (t, s)-последовательности; для *i*-й компоненты `polynomials.size()`-мерного пространства будет использован `polynomials[i-1]`-ый многочлен.

Этот конструктор следует использовать, если требуется вручную настроить однородность точек вдоль каждой оси. В отличие от второго конструктора, данный конструктор требует явного задания многочленов.

[^ наверх ^](#содержание)


#### Включение файлов в Вашу программу

При написании программ, использующих данный генератор, потребуются файлы из папки `include`.

В своём исходном коде пропишите строку

    #include "<your/path/to/our/files/>include/niederreiter2.hpp"

[^ наверх ^](#содержание)


#### Примеры

###### Пример 1. Простейший способ использования

Участок кода ниже генерирует и выводит первые 10 точек (t, s)-последовательности в трёхмерном пространстве.

    uint32_t dim    = 3;
    uint32_t amount = 10;
    sequences::Niederreiter<uint64_t, 64>   generator(dim);
    sequences::Point                        tmp_point;
    for (uint32_t i = 0; i < amount; ++i)
    {
        tmp_point = generator.get_point_real(i); // генерирует i-ую точку последовательности
        for (uint32_t j = 0; j < dim; ++j)
            std::cout << tmp_point[j] << '\t';
        std::cout << '\n';
    }

Вывод программы:

    0       0       0
    0.5     0.5     0.75
    0.75    0.25    0.3125
    0.25    0.75    0.5625
    0.375   0.375   0.875
    0.875   0.875   0.125
    0.625   0.125   0.6875
    0.125   0.625   0.4375
    0.1875  0.3125  0.515625
    0.6875  0.8125  0.265625

###### Пример 2. Замена старого генератора новым

Библиотека предусматривает случаи, когда разработчикам требуется использовать несколько генераторов с разными параметрами последовательно. Для этого можно воспользоваться оператором присваивания `=` — генератор, хранившийся в памяти до этого, уничтожается. Так, участок кода ниже генерирует 5 точек из 5 различных последовательностей в четырёхмерном пространстве.

    uint32_t                          amount   = 5;
    std::vector<sequences::BasicInt>  poly_deg = {1, 2, 3, 4};
    sequences::Niederreiter<uint32_t, 32>  generator(poly_deg);
    sequences::Point                       tmp_point;
    for (uint32_t i = 0; i < amount; ++i)
    {
        tmp_point = generator.get_point_real(i);
        for (uint32_t j = 0; j < poly_deg.size(); ++j)
            std::cout << tmp_point[j] << '\t';
        std::cout << '\n';
        poly_deg = {i+2, i+3, i+4, i+5}; // задаём новые степени многочленов
        generator = sequences::Niederreiter<uint32_t, 32>(poly_deg); // старый генератор стирается из памяти, заменяется новым
    }

Вывод программы:

    0               0               0               0
    0.75            0.875           0.9375          0.96875
    0.140625        0.0664062       0.0322266       0.0158691
    0.878906        0.938477        0.968994        0.984436
    0.0644531       0.0317383       0.0157471       0.00784302

###### Пример 3. Использование обработчиков генерации

Библиотека позволяет назначать пользовательские функции, которые будут обрабатывать каждую сгенерированную точку из последовательности. Например, участок программы ниже выводит первую координату каждой точки.

    // Каждый обработчик должен иметь такую сигнатуру
    void my_point_handler(sequences::Point const &new_point, sequences::CountInt point_i)
    {
        std::cout << "point[" << point_i << "]'s first component is " << new_point[0] << '\n';
        return;
    }
    
    <...>
    
    uint32_t                          start_i  = 3;
    uint32_t                          amount   = 5;
    std::vector<sequences::BasicInt>  poly_deg = {3, 2, 8, 5};
    sequences::Niederreiter<uint32_t, 32>  generator(poly_deg);
    sequences::Point                       tmp_point;
    generator.for_each_point_real(my_point_handler, amount, start_i); // Генерирует amount точек, начиная со start_i-ой, и посылает их в my_point_handler

Вывод программы:

    point[3]'s first component is 0.765625
    point[4]'s first component is 0.28125
    point[5]'s first component is 0.65625
    point[6]'s first component is 0.421875
    point[7]'s first component is 0.546875

[^ наверх ^](#содержание)




## Тестирование

Информация о тестировании может быть найдена в документе [TsTests Usage and Development Guide](https://github.com/jointpoints/tms-nets/blob/master/tests/README.md) (на английском языке).

[^ наверх ^](#содержание)




## Для дальнейшей разработки


#### Документация

Документация к стабильной версии программы доступна [здесь](https://jointpoints.github.io/tms-nets/) (частично на английском языке),
генерируется автоматически с помощью утилиты Doxygen при каждом изменении в ветви **master**. Для полного понимания номенклатуры и алгоритма рекомендуется ознакомиться с литературой из списка ниже.

[^ наверх ^](#содержание)


#### Рекомендуемая литература

  1. Harald Niederreiter "Low-Discrepancy and Low-Dispersion Sequences"
  2. Harald Niederreiter "Random Number Generation and Quasi-Monte Carlo Methods" (особенно глава 4)
  3. Paul Brately, Bennett L. Fox, Harald Niederreiter "Implementation and Tests of Low-Discrepancy Sequences"

[^ наверх ^](#содержание)
