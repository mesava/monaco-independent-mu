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

## Обновление M1 quality workflow (20.11.2026)

Реализованы:

- безопасное чтение headerless implicit-VR Monaco DICOM;
- локальный seed snap вместо глобального 3D EDT;
- одна загрузка CT для трёх пороговых вариантов;
- Dice, добавленные/удалённые voxels и относительный объём;
- posterior survey для strictly axial HFS;
- отдельная команда `patient-qc-png` для просмотра CT slices;
- вывод локального QA JSON.

В этом обновлении **не утверждается**, что независимый визуальный контроль
стола/подложки завершён. Новые patient-specific числа и окончательный
`.egsphant` не опубликованы до воспроизводимого прогона и проверки.

Подробнее: `docs/m1_mask_quality.md`.

## Повторный расчёт raw CT — 20.11.2026

Исходный архив с 187 CT instances повторно обработан в защищённой локальной
вычислительной среде **непосредственно из пиксельных данных DICOM**. Для
проверки чисел использован отдельный простой implicit-VR reader и
segmentation algorithm, соответствующий текущему research path:
3-D binary closing=1, 6-connected component selected by isocenter,
nearest-tissue seed-snap <=50 mm, per-slice binary hole filling.

Контроль совпадения с ранее записанным baseline при -500 HU:

- 6 435 424 voxels;
- 39 278.7109375 cm³;
- 5 898 483 voxels до per-slice hole filling;
- HU in isocenter = -821; seed snap = 2.209708691 mm.

### Три порога HU

| Показатель | -600 HU | -500 HU | -400 HU |
|---|---:|---:|---:|
| Воксели в маске | 6 641 182 | 6 435 424 | 6 181 312 |
| Объём, см³ | 40 534.5581 | 39 278.7109 | 37 727.7344 |
| Dice относительно -500 HU | 0.984265 | 1.000000 | 0.979859 |
| Добавлено вокселей к -500 HU | 205 758 | 0 | 0 |
| Исключено вокселей относительно -500 HU | 0 | 0 | 254 112 |
| HU < -1000, вокселей внутри маски | 155 | 42 | 0 |
| HU > +2009, вокселей внутри маски | 13 | 13 | 13 |
| Плотные воксели позади всей маски | 915 | 750 | 3 402 |
| Срезы с такими вокселями | 42 | 33 | 51 |

Воксели outside scanner calibration при -500 HU:
42 ниже нижней границы и 13 выше верхней, **55 всего**,
или ~0.000855% маски. Это не снимает необходимости физического анализа
их локализации (например, металл и артефакты) и чувствительности
вычисленной дозы к policy endpoint clipping.

### Критическая находка: в маску включена опора

На локальной QA-монтажной сетке CT с наложенным research contour видно, что
выбранная 3-D связная компонента включает **горизонтальные части CT-стола /
фиксационного устройства и нижнюю дугообразную конструкцию**. Эти элементы
обведены той же границей, что и тело пациента, на разных уровнях z.

**39 278.71 см³ НЕ является объёмом только пациента.**

Предыдущий `posterior_survey` проверял материал *вне* маски за её самой
задней границей. Он не способен выявить конструкцию, уже попавшую в маску
по пути 3-D connectivity. Поэтому проект теперь дополнительно использует
`survey_posterior_mask_islands`, анализирующий отдельные posterior 2-D
островки **внутри** маски, и приоритетно показывает подозрительные срезы
на локальном QA PNG.

Наличие таких 2-D островков — признак для ручной проверки, но не
автоматическая классификация стола. Отсутствие островков также не
подтверждает, что соединённая с телом опора исключена.

**Статус: M1 BODY/SUPPORT SEGMENTATION FAILED; RESEARCH_MASK_UNVALIDATED.**
Анализ представляет собой реальную оценку геометрии и HU, **не**
commissioned patient model, не расчёт дозы и не основание выбирать
clinical RED policy.

### Что делать перед egsphant

Предпочтительно получить независимо просмотренный External/BODY contour
на той же CT и в том же FrameOfReferenceUID, при необходимости
отдельную Couch/Support structure. Альтернативный путь — разработать
и отдельно валидировать segmentation patient/support на серии CT
с trusted reference contours.

Вариант «вырезать всё ниже фиксированной Y-координаты» не считается
универсально безопасным для breast board, рук и контактирующей с ними
опоры. До независимой проверки corrected mask не объявляется физически
подтверждённой.

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
