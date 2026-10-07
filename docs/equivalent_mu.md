# Per-beam multi-ROI equivalent MU

## Что именно проверяется

Для IMRT/VMAT нельзя честно утверждать, что MU плана независимо выводятся из
prescription без воспроизведения inverse optimization.

Поэтому используем термин **equivalent MU consistency check**.

Для каждого beam и ROI:

```text
MU_equiv =
    MU_TPS * D_TPS,beam(ROI) / D_MC,beam(ROI at MU_TPS)
```

и:

```text
DeltaMU[%] = 100 * (MU_equiv - MU_TPS) / MU_TPS
```

Это проверяет независимую dose-per-MU данного delivery.

## Три ROI

Проект предусматривает три малых ROI радиусом 0.25 cm:

1. reference / near-isocenter;
2. high-dose low-gradient PTV;
3. off-axis stress point.

Их конкретные координаты являются patient-specific и выбираются только после
анализа RTSTRUCT/RTDOSE.

## Почему не берём просто raw voxels

Радиус 2.5 mm может быть сравним с dose-grid spacing. Поэтому среднее считается
в непрерывном пространстве:

1. симметричное supersampling sphere;
2. trilinear interpolation dose grid;
3. среднее по samples внутри sphere.

Таким образом результат меньше зависит от случайного положения центра ROI
относительно voxel lattice.

## ROI spread

Отдельно считается:

```text
spread[% of TPS MU] =
    100 * (max(MU_equiv) - min(MU_equiv)) / MU_TPS
```

То есть spread всегда является процентом, а не разностью MU ошибочно
подписанной знаком `%`.

## Tolerances

Infrastructure принимает tolerance как конфигурацию.

Исходно для проекта зафиксированы рабочие цели:

- простые планы: ±3%;
- IMRT/VMAT: ±5%.

Однако автоматические границы YELLOW/RED не придумываются заранее. Они будут
зафиксированы только после validation dataset и локального клинического
решения.
