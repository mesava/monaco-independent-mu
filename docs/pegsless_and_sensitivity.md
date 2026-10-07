# PEGSless media и sensitivity framework

## Что реализовано

После дискретизации непрерывного Monaco Patient LUT каждый medium из
`.egsphant` получает одноимённое PEGSless-описание:

- elemental composition;
- mass fractions;
- reference bulk density;
- explicit `bremsstrahlung correction = NRC`.

Блок генерируется в формате EGSnrc:

```text
:start media definition:
    ...
    :start MEDIUM:
        elements = ...
        mass fractions = ...
        bulk density = ...
        bremsstrahlung correction = NRC
    :stop MEDIUM:
:stop media definition:
```

Порядок media **строго совпадает** с порядком medium names в `.egsphant`.

## I-value

Для дискретных промежуточных смесей I-value пока намеренно не задаётся вручную.
EGSnrc рассчитывает density-effect data on-the-fly из elemental composition,
что соответствует рекомендуемому PEGSless workflow для новых материалов в
PIRS-701.

Это позволяет не вводить неподтверждённое правило интерполяции I-value между
ICRP/NIST тканями.

## Reference density

`.egsphant` уже содержит фактическую массовую плотность каждого voxel.
PEGSless medium, однако, также требует reference bulk density.

В текущем baseline reference density определяется из reference densities двух
базовых составов линейно по material-mix fraction. Это **не подменяет** voxel
density: DOSXYZnrc использует voxel density через RHOR/RHO(MEDIUM).

Оставшаяся зависимость density-effect correction от reference density должна
быть проверена чувствительностью до commissioning.

## Sensitivity framework

Добавлены метрики для кандидатов `mixture_bins = 5, 9, 19`:

- фактическое число EGSnrc media;
- mean / p95 / max ошибки квантования material-mix fraction;
- min/max RHOR/RHO(MEDIUM);
- p95 отклонения RHOR/RHO(MEDIUM) от 1;
- доля voxels с отклонением reference-density ratio более 20%.

Эти метрики являются **структурной проверкой patient model**, а не заменой
дозовой sensitivity study.

Окончательный выбор числа media будет сделан только после запуска transport и
сравнения dose distributions для нескольких discretization schemes.

## Следующий этап

1. интегрировать PEGSless media block в DOSXYZnrc input generator;
2. в M4 прогнать dose sensitivity 5 vs 9 vs 19 bins;
3. если RHOR/RHO(MEDIUM) окажется значимым источником различий, добавить
   отдельное reference-density binning с учётом ограничения 61 medium.
