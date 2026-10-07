# Validation baseline

Для реальных Monaco exports проект разделяет два артефакта:

1. **technical fingerprint** — фактическое обезличенное описание DICOM case;
2. **expected baseline** — утверждённый набор технических свойств, который
   будущая версия parser должна продолжать интерпретировать одинаково.

Для первого case baseline находится в:

    config/validation/case001.expected.yaml

Проверка выполняется командой:

    indep-mu validate-baseline <DICOM_DIR> config/validation/case001.expected.yaml

## Что проверяется

Baseline case 001 фиксирует:

- CT matrix, spacing, orientation, HFS и 120 kV;
- количество ROI и отсутствие exported external/BODY/patient ROI;
- 25 fractions;
- 8 treatment beams;
- total MU per fraction;
- 6-MV STANDARD energy resolution;
- DYNAMIC_MLC delivery class;
- общий isocenter;
- MLCX 80 × 5 mm, span 400 mm, DICOM reference plane 349 mm;
- ASYMY reference plane 470 mm;
- MU, CP count, gantry и collimator каждого beam;
- наличие только plan-level RTDOSE в Gy;
- отсутствие per-beam RTDOSE;
- unknown physical dose quantity в стандартном RTDOSE.

## Зачем отдельная команда

Unit tests с synthetic DICOM проверяют код, но не доказывают, что parser
по-прежнему одинаково трактует реальный Monaco export.

validate-baseline позволяет прогнать именно локальный raw DICOM, не
коммитируя пациентские файлы в Git.

Сравнение является expected-subset comparison: новые диагностические поля
могут появляться без поломки baseline, а изменение уже утверждённого поля
даёт FAIL с точным путём, например:

    rtplan.beams[3].control_points:
        expected=29, actual=30

case_id и текстовый status являются metadata baseline-файла и не
сравниваются с DICOM.

## Что baseline не доказывает

PASS подтверждает только стабильность DICOM interpretation.

Он не подтверждает:

- корректность CT-derived external mask;
- physical Agility head geometry;
- absolute MC calibration;
- dose-to-medium/dose-to-water semantics;
- точность independent MC dose.

Эти пункты имеют отдельные validation gates.
