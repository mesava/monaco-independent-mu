# RTDOSE geometry

RT Dose загружается как физическая 3-D сетка в DICOM patient coordinate system.

## DoseGridScaling

Raw Pixel Data умножается на `DoseGridScaling`. Единицы сохраняются из
`DoseUnits`; для будущего абсолютного Monaco-vs-MC сравнения ожидается `GY`.

## GridFrameOffsetVector

Реализованы оба варианта DICOM PS3.3 C.8.8.3.2.

### Relative — основной современный вариант

Если первый элемент GridFrameOffsetVector равен 0:

```text
P_frame(k) =
    ImagePositionPatient
    + GridFrameOffsetVector[k] * normal

normal = row_direction × column_direction
```

Такой вариант поддерживает произвольную ориентацию dose grid.

### Legacy absolute patient-z

Если orientation строго `(1,0,0,0,1,0)`, а первый элемент vector равен
третьему компоненту ImagePositionPatient, значения трактуются как абсолютные
patient-z координаты.

Для обlique grid этот legacy режим запрещён.

## Безопасность

Проверяются:

- NumberOfFrames vs GridFrameOffsetVector length;
- строгая монотонность offsets;
- pixel array shape;
- обязательные DoseUnits / DoseType / DoseSummationType;
- FrameOfReferenceUID.

Dose geometry не предполагается совпадающей с CT grid. Для 3D comparison
регистрация и интерполяция выполняются отдельным модулем в patient coordinates.
