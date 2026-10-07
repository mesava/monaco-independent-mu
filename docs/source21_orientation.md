# DICOM IEC → DOSXYZnrc Source 21 orientation

## Реализованный mapper

В проект добавлен `ZhanHfsSource21OrientationMapper`.

Он реализует компактное преобразование DICOM gantry/couch/collimator →
DOSXYZnrc `theta/phi/phicol`, опубликованное Zhan, Jiang & Osei и
впоследствии опубликованное автором Lixin Zhan как reference code.

Для HFS:

```text
gamma = gantry
rho   = couch / PatientSupportAngle
col   = collimator

theta = acos(-sin(gamma) sin(rho))

phi = atan2(
    -cos(gamma),
     sin(gamma) cos(rho)
)

couch_projection =
    atan2(
        -sin(rho) cos(gamma),
         cos(rho)
    )

phicol =
    pi - [
        (col - pi/2)
        + couch_projection
    ]
```

После этого `phi` и `phicol` приводятся к 0…360°.

## Независимые regression checks

При couch=0 mapper обязан сводиться к ранее зафиксированной coplanar HFS
формуле:

```text
theta  = 90°
phi    = -90° + gantry
phicol = -90° - collimator
```

Это проверяется unit tests для cardinal и произвольных coplanar angles.

Дополнительно в тестах зафиксированы несколько non-coplanar vectors,
независимо рассчитанных по опубликованной формуле.

## Singularities

Reference implementation автора слегка изменяет gantry и couch
(`×0.999999`) при комбинациях:

```text
gantry = 90° или 270°
couch  = 90° или 270°
```

чтобы обойти singular spherical-angle decomposition.

Для клинического независимого калькулятора такое скрытое изменение delivery
angles неприемлемо.

Поэтому production mapper **останавливает расчёт** на этих точных комбинациях.
Если позже понадобится их поддержка, будет реализована отдельная
rotation-matrix representation без perturbation и проверена геометрически.

## Patient position

Текущая реализация разрешает только:

```text
PatientPosition = HFS
```

FFS/HFP/FFP/decubitus нельзя просто пропустить через ту же формулу, поскольку
связь IEC patient-support coordinates с DICOM patient coordinates меняется.

До отдельной валидации таких положений mapper их отвергает.

## Источники

- Zhan L, Jiang R, Osei EK. *Beam coordinate transformations from DICOM to
  DOSXYZnrc.* Phys Med Biol. 2012;57:N513-N523.
- Lixin Zhan, PYCOM reference implementation `dicom2dosxyz.py`.
- SlicerRT `EGSnrcUtil.dcm2dosxyz()`, реализация метода Zhan.
- NRC/EGSnrc 2026 course, DOSXYZnrc coordinate-transformation laboratory.
