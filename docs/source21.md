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

Добавлен `ZhanHfsSource21OrientationMapper`, реализующий опубликованное
преобразование Zhan/Jiang/Osei для HFS.

Проверки включают:

- gantry 0/90/180/270;
- collimator rotations;
- согласие с coplanar reference equations;
- несколько non-coplanar regression vectors.

Точные сочетания gantry=90/270 и couch=90/270 считаются singular и
останавливают расчёт вместо скрытого perturbation angles.

Patient positions кроме HFS пока запрещрещены до отдельной валидации.

Подробности: `docs/source21_orientation.md`.

## Следующий шаг

Остаётся геометрическая end-to-end validation source orientation уже в самом
DOSXYZnrc/BEAMnrc на контрольных планах. После этого mapper можно перевести из
research-validated в commissioned transport configuration.
