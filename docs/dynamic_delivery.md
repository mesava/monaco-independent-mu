# Dynamic delivery semantics

## Что гарантирует DICOM

Control Point описывает состояние delivery device в конкретной точке процесса.
Rotation Direction относится к сегменту **после** данного control point.
Cumulative Meterset Weight монотонно растёт.

Однако DICOM PS3.3 C.8.8.14.5 прямо указывает, что стандарт **не делает
предположений о поведении machine parameters между control points**.

Следовательно, линейная интерполяция MLC/jaws/gantry между CP не может быть
неявным свойством DICOM parser.

## Правило проекта

Parser сохраняет только то, что сообщает DICOM.

Monte Carlo transport обязан получить отдельный
`DeliveryInterpolationPolicy` с:

- идентификатором policy;
- evidence/validation statement;
- явным алгоритмом интерполяции.

Текущий реализованный алгоритм — linear in Cumulative Meterset Weight. Он
предназначен для будущего BEAMnrc sequence-file workflow и должен быть
валидирован для Elekta/Monaco delivery до клинического использования.

## Rotation Direction

Используется DICOM machine-rotation semantics:

- `NONE`: движения нет;
- `CC`: движение в сторону увеличения IEC angle;
- `CW`: движение в сторону уменьшения IEC angle;
- одинаковые start/end angles при CW/CC означают полный оборот 360°.

## Почему это важно

Такой дизайн предотвращает скрытую подмену:

```text
DICOM control points
        !=
гарантированная непрерывная траектория машины
```

и позволяет позже сравнить planned sequence с Elekta delivery logs/iCOM при
валидации динамической модели.
