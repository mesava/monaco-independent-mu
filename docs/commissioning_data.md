# Commissioning data model

## Цель

Commissioning dataset является **независимым измерительным входом** beam model.
Он отделён от Monaco TPS tuning parameters и от доз, рассчитанных самой TPS.

Поддерживаются:

- absolute reference calibration;
- output factors;
- PDD;
- in-plane / cross-plane profiles.

## Абсолютная калибровка

Независимая абсолютная калибровка задаётся только в размерно однозначном виде:

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

Поле `delivered_mu` является обязательным. Проект никогда не угадывает,
соответствует ли значение дозы 1 MU, 100 MU или другому числу MU.

## Monaco TPS reference dose

Отдельно хранится `config/machine/versa_hd/monaco_reference.yaml`.

Это **не commissioning calibration** и не используется как независимый источник
абсолютной нормировки Monte Carlo. Это regression/sanity reference для
сопоставления с Monaco в зафиксированной геометрии:

- поле 10×10 cm²;
- SSD 90 cm;
- глубина 10 cm;
- 100 MU.

Подтверждённые значения Monaco 6.1.4:

| Energy model | Monaco dose for 100 MU |
|---|---:|
| 6MV | 0.995 Gy |
| 6FFF | 0.994 Gy |
| 10MV | 1.000 Gy |
| 10FFF | 1.001 Gy |

Соответствующие TPS dose-per-MU:

- 6MV: 0.00995 Gy/MU;
- 6FFF: 0.00994 Gy/MU;
- 10MV: 0.01000 Gy/MU;
- 10FFF: 0.01001 Gy/MU.

Эти числа полезны для regression testing и сравнения, но использовать их для
нормировки независимого MC означало бы частично замкнуть независимую проверку
на TPS.

## Output factors

Output factor хранится вместе с двумя размерами поля. Для reference field
10×10 cm² предусмотрена отдельная проверка близости OF к 1.

## PDD и profiles

Сохраняется полная геометрия измерения:

- field size;
- SSD;
- depth;
- ось профиля;
- исходные coordinate/value arrays.

Нормировка кривых не должна уничтожать исходные commissioning data.
Импорт конкретных ASC/TXT форматов выполняется отдельным адаптером.

## Что ещё требуется для независимой абсолютной MC-нормировки

Нужен **измеренный** reference dose на ускорителе в однозначно заданной
референсной геометрии с известным количеством MU. Он может совпадать с
клиническим calibration point, но должен происходить из измерения, а не из
Monaco RT Dose.

Monaco LT-параметры также не являются заменой измеренным PDD/profiles/OF.
