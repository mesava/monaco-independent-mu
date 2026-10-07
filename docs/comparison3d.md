# M6 — 3D dose comparison

## Общий принцип

Reference — физический RTDOSE Monaco.

Evaluation — независимо нормированная MC dose in Gy.

Никакой best-fit normalization TPS ↔ MC внутри comparison layer нет.

## Patient-coordinate MC adapter

PatientMcDoseGridGy связывает DOSXYZnrc 3ddose, egsphant axis directions и
FrameOfReferenceUID.

Любая DICOM patient coordinate переводится в axis basis MC phantom и
семплируется trilinear interpolation. Поэтому TPS и MC ROI используют одни и
те же физические координаты.

## Voxel-wise dose difference

compare_on_rtdose_grid() оставляет Monaco RTDOSE в качестве reference grid и
семплирует только evaluation MC на центры TPS voxels.

Основная величина:

    DeltaD = D_MC - D_TPS

Relative difference вычисляется только внутри dose-threshold mask.

Рабочий baseline threshold — D_TPS не менее 10% от максимальной reference dose,
но окончательная clinical policy остаётся конфигурируемой.

## DVH

Структурные dose samples берутся на центрах CT voxels внутри RTSTRUCT mask.

Реализованы:

- D98;
- D95;
- D50;
- D2;
- Dmean;
- Dmin;
- Dmax.

Используется стандартная percentile semantics:

    D98 = 2nd percentile
    D95 = 5th percentile
    D50 = median
    D2  = 98th percentile

Позже для OAR добавятся endpoint-specific volume metrics и D0.1cc/D2cc там,
где это клинически оправдано.

## Gamma

Gamma намеренно пока не смешан с dose-difference implementation.

Перед production implementation должны быть явно зафиксированы:

- global vs local criterion;
- normalization dose;
- dose threshold;
- interpolation/search implementation;
- validation against a trusted reference implementation.

Рабочая цель проекта остаётся primary 3%/3 mm и diagnostic PTV-core 2%/2 mm.
2%/2 mm на старте не является автоматическим clinical failure criterion.
