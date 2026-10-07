# Commissioning data model

## Цель

Commissioning dataset является независимым измерительным входом beam model.
Он отделён от Monaco TPS tuning parameters.

Поддерживаются:

- absolute reference calibration;
- output factors;
- PDD;
- in-plane / cross-plane profiles.

## Абсолютная калибровка

Калибровка задаётся только в размерно однозначном виде:

```yaml
field: 100 x 100 mm
ssd_mm: ...
depth_mm: ...
delivered_mu: ...
measured_dose_gy: ...
```

После этого:

```text
dose_per_MU = measured_dose_gy / delivered_mu
```

Проект **никогда не предполагает автоматически**, что измерение 0.995 Gy
соответствует 100 MU. Даже если это наиболее вероятная клиническая геометрия,
число MU должно быть явно подтверждено в исходных данных.

Это защищает от ошибки в 100 раз при абсолютной нормировке MC.

## Output factors

Output factor хранится вместе с двумя размерами поля. Для reference field
10x10 cm² предусмотрена отдельная проверка близости OF к 1.

## PDD и profiles

Сохраняется полная геометрия измерения:

- field size;
- SSD;
- depth;
- ось профиля;
- исходные coordinate/value arrays.

Нормировка кривых не должна уничтожать исходные commissioning data.
Импорт конкретных ASC/TXT форматов будет отдельным адаптером.

## Что пока не зафиксировано

Предоставленные ранее reference-dose значения для 6MV / 6FFF / 10MV / 10FFF
не внесены в рабочий commissioning YAML, потому что число доставленных MU
должно быть подтверждено явно.

Monaco LT-параметры также не являются заменой измеренным PDD/profiles/OF.
