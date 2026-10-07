# Technical DICOM fingerprint

Команда:

    indep-mu fingerprint <DICOM_DIR>

формирует технический JSON fingerprint одного DICOM case.

По умолчанию результат **не содержит**:

- PatientName;
- PatientID;
- raw DICOM UIDs;
- plan/beam names;
- ROI names.

UID сохраняется только как короткий SHA-256 fingerprint для проверки
внутренней согласованности без копирования исходного идентификатора.

ROI names можно включить отдельным флагом:

    --include-roi-names

Это сделано явно, потому что произвольное название структуры теоретически
может содержать персональный текст.

Fingerprint включает:

- CT matrix/spacing/orientation/KVP;
- observed HU range и whole-CT counts вне scanner calibration;
- plan/fraction/beam/MU/control-point metadata;
- delivery class;
- Agility delivery geometry;
- plan/transport preflight;
- RTDOSE geometry, units/type/summation и dose range.

## Validation case 001

Проверенный metadata baseline хранится в:

    config/validation/case001.expected.yaml

Он не содержит исходных DICOM UIDs или персональных данных. Цель файла —
регрессионно проверить, что будущие изменения parser не меняют уже
подтверждённую интерпретацию первого реального Monaco export.

Pixel-level HU/body statistics в expected file пока намеренно отсутствуют:
они будут добавлены только после воспроизводимого запуска текущего pipeline на
raw DICOM в рабочей вычислительной среде.
