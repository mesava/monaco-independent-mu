# Versa HD / Agility Y-jaw geometry: current evidence and blocker

## Что подтверждает реальный DICOM

Validation case 001 экспортирует для ASYMY:

    SourceToBeamLimitingDeviceDistance = 470 mm

DICOM определяет (300A,00BA) как расстояние от radiation source до
beam-limiting device. Стандарт не определяет этим атрибутом физические
upstream/downstream surfaces jaw block.

Следовательно, значение 470 mm полезно как delivery reference plane, но само
по себе не задаёт zmin/zmax, требуемые SYNCJAWS.

## Почему не используем Monaco LT как физическую геометрию автоматически

В Monaco reference parameters встречается TJawPlanePosition около 432 mm.
Это tuned TPS/model parameter. Без документации о его точной физической
семантике нельзя объявлять его upstream face jaw и выводить толщину как
2 × (470−432).

Такая арифметика выглядит правдоподобно, но была бы скрытым
manufacturer-model assumption.

## Что говорит открытая литература

Опубликованный Versa HD EGSnrc model использует SYNCJAWS для Y jaws, но
авторы прямо указывают, что подробные dimensions treatment head были получены
от Elekta и являются proprietary.

Публикация Infinity/Agility аналогично строит детальный head model на основе
manufacturer data и использует synchronized lower jaws.

Таким образом открытая литература подтверждает методику, но не даёт
достаточной публичной геометрии для однозначного zmin/zmax нашей Y-jaw.

## Безопасная стратегия проекта

До появления однозначной geometry evidence:

1. DICOM 470 mm хранится как reference-plane observation;
2. SYNCJAWS zmin/zmax остаются blocked;
3. Monaco LT parameters не конвертируются в physical dimensions;
4. patient transport не запускается с выдуманной jaw thickness.

Допустимые способы закрыть blocker:

- Elekta physical head geometry / engineering documentation;
- уже валидированная независимая Versa HD BEAMnrc geometry;
- measurement-driven effective jaw model, если его параметры оказываются
  идентифицируемыми по профилям/transmission и проходят независимую
  validation series.

Последний вариант должен называться effective model, а не manufacturer geometry.

## Источники

- DICOM PS3.3 RT Beams Module, tag (300A,00BA).
- Paschal et al. JACMP 2022, DOI 10.1002/acm2.13715.
- Gholampourkashi et al. JACMP 2019, DOI 10.1002/acm2.12485.
