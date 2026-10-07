# EGSnrc: voxel density scaling и density-effect caveat

## Что делает DOSXYZnrc с density из .egsphant

Для каждого geometry region EGSnrc допускает плотность `RHOR`, отличную от
номинальной плотности среды `RHO(MED)`.

Официальная документация EGSnrc описывает масштаб:

```text
RHOF = RHOR(IRL) / RHO(MEDIUM)
```

и указывает, что cross sections в регионе соответствующим образом
масштабируются.

Источник:
https://github.com/nrc-cnrc/EGSnrc/blob/master/HEN_HOUSE/doc/src/pirs701-egsnrc/inputs/new_um_app2.tex

В стандартных macros это реализовано как:

```text
$SET-RHOF → RHOF=RHOR(IRL)/RHO(MEDIUM)
```

Источник:
https://github.com/nrc-cnrc/EGSnrc/blob/master/HEN_HOUSE/src/egsnrc.macros

Следовательно, voxel-wise mass density из нашего PatientModel действительно
может записываться в `.egsphant`, даже если несколько voxel используют одну
и ту же дискретную elemental composition.

## Ограничение

EGSnrc отдельно предупреждает, что такое density scaling **не является
абсолютно точным**, потому что density-effect correction в electron stopping
powers не масштабируется вместе с RHOR.

Для высокоточной задачи документация рекомендует при необходимости определять
несколько media одного состава с разными nominal densities.

Это принципиально важно для нашего проекта.

## Решение в архитектуре

Мы разделяем две независимые дискретизации:

1. **composition discretization**
   - аппроксимирует непрерывный Monaco MaterialMix;
   - текущие sensitivity candidates: 5 / 9 / 19 mixture intervals.

2. **density-reference discretization**
   - несколько EGSnrc media могут иметь одинаковый elemental composition,
     но различные nominal/local densities;
   - количество density sub-bins пока НЕ выбрано.

Параметр density-reference discretization будет выбран только после MC
sensitivity study. До этого pegsless media generator не считается
commissioning-final.

## Что проверять

Минимальная sensitivity matrix:

- composition bins: 5, 9, 19;
- density-reference variants:
  - one nominal density per composition;
  - additional density bins в наиболее широких диапазонах;
- homogeneous soft tissue;
- low-density lung;
- cortical bone;
- heterogeneous slab;
- Patient Validation Case 001.

Основная метрика: изменение absolute dose и 3D dose при увеличении
дискретизации. Если результат стабилизируется раньше технического максимума,
используется минимальная стабильная конфигурация.

Таким образом, лимит 61 medium legacy egsphant рассматривается как
техническое ограничение backend, а не как физический критерий.
