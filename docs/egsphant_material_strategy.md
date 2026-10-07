# Дискретизация Monaco MaterialMix для EGSnrc / DOSXYZnrc

## Почему дискретизация вообще нужна

Monaco Patient LUT задаёт непрерывное смешивание соседних материалов по
mass density. Например, между 1.08 и 1.85 g/cm³ используется смесь
MuscleSkeletalIcrp и BoneCorticalIcrp.

DOSXYZnrc в CT-режиме хранит для каждого voxel:

1. **номер среды**;
2. **плотность voxel**.

То есть плотность может меняться voxel-by-voxel, но элементный состав среды
остаётся дискретным. Поэтому непрерывный Monaco MaterialMix нельзя буквально
записать в legacy `.egsphant` без конечного набора промежуточных сред.

Официальный DOSXYZnrc manual:
https://nrc-cnrc.github.io/EGSnrc/doc/pirs794-dosxyznrc.pdf

## Технический предел legacy .egsphant

В стандартном `ctcreate.mortran` используется строка кодирования длиной 62
символа и вывод `encoding(medium_index + 1)`. Поэтому стандартный формат
имеет максимум **61 используемый номер среды**.

Исходный код:
https://github.com/nrc-cnrc/EGSnrc/blob/master/HEN_HOUSE/omega/progs/ctcreate/ctcreate.mortran

Кроме того, в DOSXYZnrc `$MXMED` по умолчанию равен 5; для нашей задачи
reference build должен быть перекомпилирован с увеличенным значением.

## Что это означает для нашего Patient LUT

В полном Patient LUT есть 7 базовых чистых материалов:

- DryAir;
- MuscleSkeletalIcrp;
- AdiposeTissueICRP;
- BoneCorticalIcrp;
- Titanium;
- StainlessSteel316;
- Lead.

Непрерывные смеси нужны только для трёх уникальных пар состава:

1. DryAir ↔ MuscleSkeletalIcrp;
2. MuscleSkeletalIcrp ↔ AdiposeTissueICRP;
3. MuscleSkeletalIcrp ↔ BoneCorticalIcrp.

Переходы Bone→Titanium, Titanium→Steel и Steel→Lead происходят при одинаковых
threshold density и являются резкими, а не интерполированными диапазонами.

Если использовать **19 интервалов смеси**, каждая пара имеет 18 внутренних
композиций:

```text
7 pure + 3 × 18 mixed = 61 media
```

Это максимальная сетка смеси, которая ещё помещается в стандартное
62-символьное кодирование `.egsphant`.

## Почему 19 пока НЕ является клиническим параметром

Выбор количества bins влияет на элементный состав voxel и потенциально на дозу.
Поэтому параметр не фиксируется по удобству формата.

Сначала сравниваются как минимум:

- 5 bins;
- 9 bins;
- 19 bins.

На тестах:

1. water/soft-tissue;
2. lung-like low density;
3. adipose↔muscle boundary;
4. muscle↔cortical bone;
5. пациентский кейс.

Критерий выбора — отсутствие практически значимого изменения дозы при
дальнейшем увеличении разрешения material mixture.

## Базовые составы

Для текущих тканевых материалов используются независимые reference
compositions NIST STAR с ICRP-наименованиями:

- Air, Dry;
- Adipose Tissue (ICRP);
- Muscle, Skeletal (ICRP);
- Bone, Cortical (ICRP).

Это не экспорт внутреннего material database Monaco. Названия LUT Monaco
сопоставляются с независимыми опубликованными ICRP/NIST compositions.

NIST STAR:
https://physics.nist.gov/PhysRefData/Star/Text/table2.html

EGSnrc также содержит стандартные AIR и ICRP cortical bone definitions в
`HEN_HOUSE/pegs4/data/material.dat`.

## Статус

Сейчас реализованы:

- непрерывное Monaco material assignment;
- компактное хранение lower/upper material + mix fraction;
- канонизация обратных пар (Muscle→Adipose и Adipose→Muscle);
- конфигурируемая дискретизация;
- контроль лимита 61 medium;
- reference elemental compositions для четырёх тканевых сред.

Следующий этап — генератор media definition + `.egsphant` writer и
sensitivity tests 5/9/19 bins.
