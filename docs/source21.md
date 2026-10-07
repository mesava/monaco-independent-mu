# DOSXYZnrc source 21: synchronization layer

## Что реализовано

Создан отдельный слой для DOSXYZnrc `ISOURCE=21`.

Он формирует записи control points:

```text
xiso yiso ziso theta phi phicol dsource muIndex
```

и проверяет, что `muIndex` **точно совпадает** с synchronized sequences
`SYNCMLCE` и `SYNCJAWS`.

Это важно: BEAMnrc/DOSXYZnrc использует один и тот же fractional monitor-unit
index для одновременного движения source orientation, MLC и jaws.

## Isocenter

DICOM `IsocenterPosition` задан в patient coordinates (mm).

Перед передачей DOSXYZnrc точка проецируется на три axis vectors конкретного
`.egsphant`. Поэтому преобразование не предполагает автоматически axial HFS
orientation.

## dsource

`dsource` **не выводится из RTPLAN автоматически**.

Он зависит от того, где находится origin/scoring plane используемой BEAMnrc
transport geometry. Поэтому значение должно поступать из валидированной MC
head configuration.

## theta / phi / phicol

Автоматическое преобразование IEC/DICOM angles пока намеренно не включено.

Причина: DICOM gantry/couch/collimator angles и DOSXYZnrc
`theta/phi/phicol` имеют разные coordinate conventions; опубликованные
работы специально посвящены этому преобразованию.

Код требует объект `Source21OrientationMapper`. Пока не выбран и не
верифицирован конкретный mapper, production source-21 input создать нельзя.

Это сделано специально, чтобы формула для couch/gantry не оказалась скрытым
непроверенным предположением.

## Следующий шаг

Orientation mapper будет валидироваться минимум на cardinal geometries:

- gantry 0/90/180/270;
- collimator 0/90/270;
- затем non-zero couch;
- визуальная/геометрическая проверка field direction в DOSXYZnrc;
- round-trip against known test plans.

После этого можно собирать полный per-beam DOSXYZnrc input.
