# Physical head geometry readiness

## Зачем отдельная конфигурация

DICOM RTPLAN задаёт delivery geometry, но не является полной физической
геометрией treatment head.

Monaco LT-параметры также не переносятся автоматически в BEAMnrc: tuned TPS
parameter не становится физическим размером только потому, что имеет похожее
название.

Поэтому physical head model хранится отдельно:

`config/machine/versa_hd/head_geometry.example.yaml`.

## Что пока намеренно null

Для Agility rounded-tip model нужны как минимум:

- MLC zmin/zmax;
- leaf-tip radius;
- cylinder-axis position (CIL);
- bank tilt (LBROT);
- DICOM bank identity.

Для jaws нужны физические zmin/zmax и bank identity.

Эти поля оставлены `null`, пока не появится источник/измерение с однозначной
физической интерпретацией.

## Два уровня готовности

### Jaw geometry

Для source-focused jaws численно полной конфигурации достаточно, чтобы
сформировать SYNCJAWS front/back coordinates. Затем она всё равно проходит
commissioning validation.

### Agility rounded tip

Даже наличие всех чисел **не снимает блокировку**.

SYNCMLCE ENDTYPE=0 требует координату origin цилиндра rounded leaf tip, тогда
как DICOM предоставляет delivery leaf position в IEC beam-limiting-device
coordinates. Нужен отдельно проверенный transform между этими величинами.

До его валидации `ready_for_rounded_agility_syncmlce = false`.

Это сознательный safety barrier: проект не должен выдавать синтаксически
правильный sequence file с физически неподтверждённой геометрией.
