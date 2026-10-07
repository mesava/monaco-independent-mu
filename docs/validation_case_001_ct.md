# Validation case 001 — CT / RTSTRUCT patient model

> **Статус: REAL DICOM METADATA VERIFIED; EXTERNAL MASK RESEARCH-ONLY.**

Полный Monaco DICOM export повторно предоставлен 2026-10-07 и техническая
интерпретация CT/RTSTRUCT/RTPLAN/RTDOSE зафиксирована в
`docs/validation_case_001_dicom.md`.

Этот файл описывает именно patient-model gate M1.

## Проверенный CT baseline

Из raw DICOM подтверждено:

- 187 CT slices;
- matrix 512 × 512;
- PixelSpacing 1.5625 × 1.5625 mm;
- slice spacing 2.5 mm;
- axial HFS;
- KVP 120 kV;
- RescaleSlope 1;
- RescaleIntercept -1024 HU;
- observed whole-CT HU range -3024 ... +3071 HU.

Scanner-specific calibration, используемая проектом:

    DICOM3.DRT120kV
    measured HU range = -1000 ... +2009 HU

## RTSTRUCT

В RTSTRUCT 14 клинических ROI, но нет отдельного:

- External;
- BODY;
- patient.

Поэтому production-oriented path с trusted external contour для этого case
недоступен.

Нельзя использовать PTV/грудь как внешнюю маску и нельзя молча считать весь
CT FOV пациентом.

## CT-derived research external mask

Для исследовательского продолжения M1 реализован отдельный явный path:

1. threshold CT;
2. 3-D binary closing;
3. connected-component selection;
4. общий treatment isocenter как seed;
5. если seed находится в low-density cavity — nearest tissue snap с maximum
   distance;
6. per-slice enclosed-hole filling.

Для первого реального прогона:

    threshold = -500 HU
    closing iterations = 1
    isocenter = (-71.3, -84.5, +100.0) mm

Получено:

- HU в seed voxel: -821 HU;
- seed находится внутри low-density lung region;
- nearest-tissue snap distance: 2.210 mm;
- selected component before fill: 5,898,483 voxels;
- final mask: 6,435,424 voxels;
- volume: 39,278.71 cm3;
- z extent: -210 ... +250 mm;
- mask не касается CT volume border;
- HU inside filled mask: -1024 ... +3071 HU.

Статус этой маски:

    RESEARCH_MASK_DERIVED_NOT_VALIDATED

## Почему M1 ещё не закрыт полностью

Результат геометрически правдоподобен, но необходимо отдельно исключить
contamination CT couch/support.

Также mask содержит HU вне измеренной DRT120kV calibration table с обеих
сторон. Проект поэтому требует:

- явный подсчёт HU < -1000;
- явный подсчёт HU > +2009;
- явную out-of-range policy;
- никакого скрытого clipping.

Для этого добавлена команда:

    indep-mu patient-diagnostics <DICOM_DIR>

По умолчанию она сравнивает -600/-500/-400 HU и показывает volume,
seed-snap distance, out-of-calibration counts и border contact.

## Gate перед первым patient .egsphant

Перед тем как считать .egsphant Case 001 validation artifact, требуется:

1. threshold sensitivity;
2. posterior support/couch contamination check;
3. out-of-calibration counts;
4. решение по endpoint clipping;
5. sensitivity выбора mixture_bins;
6. только после этого generate/hash .egsphant + PEGSless media.

Таким образом CT parsing / HU calibration / patient-model код реализован,
но реальный M1 считается закрытым только после валидации external-mask и
material discretisation.
