# Legacy .egsphant writer

Проект содержит собственный writer/read-back validator стандартного текстового
формата `.egsphant`, используемого DOSXYZnrc.

## Encoding medium numbers

Используется точная строка из официального
`dosxyznrc_user_macros.mortran`:

```text
0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz
```

DOSXYZnrc восстанавливает номер среды как `index(encoding, character)-1`.
Поэтому:

- medium 1 → `1`;
- medium 9 → `9`;
- medium 10 → `A`;
- medium 61 → `z`.

Источник:
https://github.com/nrc-cnrc/EGSnrc/blob/master/HEN_HOUSE/user_codes/dosxyznrc/dosxyznrc_user_macros.mortran

## Координаты

Legacy `.egsphant` не хранит DICOM orientation matrix. Поэтому writer строит
локальную правую систему координат:

- x — вдоль первого вектора ImageOrientationPatient, т.е. увеличения column;
- y — вдоль второго вектора, т.е. увеличения row;
- z — cross(x,y), нормаль к CT slices.

Voxel boundaries записываются в сантиметрах, как требует DOSXYZnrc.

Связь локальной системы с DICOM patient frame должна храниться отдельно и
использоваться будущим beam-coordinate transform. Потеря orientation при
экспорте в `.egsphant` не допускается на уровне workflow metadata.

## Round-trip

Unit test выполняет:

```text
synthetic CT + PatientModel
        ↓
material discretization
        ↓
write .egsphant
        ↓
read .egsphant
        ↓
compare media map / density / boundaries
```

Это проверяет сериализацию, но ещё не является Monte Carlo validation.
