# Coplanar HFS orientation reference

До реализации общей матричной трансформации DICOM IEC → DOSXYZnrc в проект
добавлен **только reference helper**, не production mapper.

Для классического coplanar HFS случая в литературе используется:

```text
theta  = 90°
phi    = -90° + gantry
phicol = -90° - collimator
```

с приведением углов к диапазону 0…360°.

Например при collimator=0:

| DICOM gantry | DOSXYZ theta | DOSXYZ phi | DOSXYZ phicol |
|---:|---:|---:|---:|
| 0° | 90° | 270° | 270° |
| 90° | 90° | 0° | 270° |
| 180° | 90° | 90° | 270° |
| 270° | 90° | 180° | 270° |

Helper жёстко отвергает non-zero couch.

## Почему это пока не Source21OrientationMapper

Для non-coplanar beams couch rotation меняет не только azimuth, но и
beam-axis orientation; `phicol` также требует корректного projection couch
rotation в collimator plane. Публикации отдельно отмечают, что эта
трансформация нетривиальна.

Поэтому простой coplanar набор формул используется только как набор
регрессионных контрольных точек для будущей общей rotation-matrix реализации.

Источники:

- Zhan L, Jiang R, Osei EK. *Beam coordinate transformations from DICOM to
  DOSXYZnrc.* Phys Med Biol. 2012;57:N513-N523.
- Ottosson R. *Monte Carlo Treatment Planning for Advanced Radiotherapy*,
  section on the MCTP workflow and coplanar transform.
- NRC DOSXYZnrc user manual, source 21 definitions.
