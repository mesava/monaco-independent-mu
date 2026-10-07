# DICOM RTDOSE parser

## Назначение

RTDOSE используется как TPS reference dose для 3D comparison и, при наличии
per-beam dose export, как numerator для equivalent-MU check.

Он **не является входом независимого Monte Carlo транспорта**.

## Dose scaling и единицы

Pixel Data преобразуется строго по DICOM:

```text
dose = stored_pixel_value × DoseGridScaling
```

`DoseUnits` сохраняется из DICOM. Для абсолютного сравнения и ROI-dose API
требуется:

```text
DoseUnits = GY
```

`RELATIVE` никогда автоматически не интерпретируется как Gy.

## Пространственная геометрия

Сохраняются и проверяются:

- ImagePositionPatient;
- ImageOrientationPatient;
- PixelSpacing;
- GridFrameOffsetVector;
- FrameOfReferenceUID.

Внутренний массив имеет порядок:

```text
dose[frame, row, column]
```

Локальные оси dose grid:

- x: направление увеличения DICOM column index;
- y: направление увеличения row index;
- z: cross(x, y).

Функция `patient_to_local_mm()` преобразует произвольную DICOM patient
coordinate в эту basis. Поэтому ROI comparison не зависит от предположения,
что RTDOSE обязательно axial/HFS.

## GridFrameOffsetVector

Реализованы обе формы DICOM PS3.3 C.8.8.3.2.

**Relative option:** первый элемент равен 0; offsets откладываются от
ImagePositionPatient вдоль normal dose grid.

**Legacy absolute patient-z option:** допускается только для
`ImageOrientationPatient = (1,0,0,0,1,0)`, когда первый GridFrameOffsetVector
совпадает с z координатой ImagePositionPatient.

Неоднозначное non-zero значение для oblique grid является ошибкой, а не
поводом для эвристики.

## Continuous dose sampling

`sample_patient_points_gy()`:

1. принимает точки в DICOM patient coordinates, mm;
2. переводит их в локальные x/y/z dose-grid coordinates;
3. выполняет trilinear interpolation;
4. возвращает Gy.

Это является основой для spherical ROI radius 0.25 cm и позднее для
resampling/gamma, но interpolation не заменяет проверку FrameOfReferenceUID.

## RTPLAN / beam references

Parser сохраняет:

- Referenced RT Plan SOP Instance UID;
- Referenced Beam Number, когда ссылка присутствует в RTDOSE.

Перед per-beam MU-check RTDOSE должен быть однозначно связан с соответствующим
beam; имя файла или порядок файлов для этого не используются.
