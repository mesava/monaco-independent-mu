# DICOM preflight

Перед Monte Carlo каждый RTPLAN проходит отдельную проверку поддержки.

## Автоматическая классификация delivery

По изменениям control-point states beam классифицируется как:

- `STATIC`;
- `DYNAMIC_JAW`;
- `DYNAMIC_MLC`;
- `ARC`;
- `VMAT`.

Для определения динамики углов используется circular difference, поэтому
переход через 0/360 градусов не создаёт ложного большого изменения.

## Текущие hard errors

До реализации соответствующей физики расчёт блокируется при:

- RadiationType != PHOTON;
- wedge / compensator / bolus / block;
- energy switching внутри одного beam;
- отсутствующем или меняющемся IsocenterPosition;
- dynamic couch;
- dynamic collimator;
- BeamMeterset <= 0.

Это deliberately conservative policy: неподдерживаемая геометрия не должна
молча рассчитываться упрощённо.

## Warnings

Пока предупреждением, а не ошибкой являются:

- TreatmentDeliveryType != TREATMENT;
- отсутствие MLC.

Preflight отделён от dose QA. Ошибка preflight означает `CALCULATION INVALID`,
а не клинический RED из-за расхождения доз.
