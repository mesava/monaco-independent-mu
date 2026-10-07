# M3 — Beam model strategy: Versa HD / Agility

## Что считаем независимой моделью

Machine commissioning measurements можно использовать как вход независимой
модели. Внутренние tuned parameters Monaco допустимы как reference/audit data,
но не должны автоматически становиться физической геометрией Monte Carlo.

Например, параметры вида `LeafGap`, `LeafOffset`,
`InterLeafLeakage` сохраняются отдельно как **Monaco reference parameters**.
Без документации нельзя объявлять их буквальной шириной зазора, вероятностью
утечки или полной геометрией Agility.

## DICOM delivery geometry

Из RTPLAN независимо получаем:

- 80 leaf pairs для Agility — проверяется, а не hard-code'ится в parser;
- LeafPositionBoundaries;
- MLCX/MLCY device identity;
- jaws;
- source-to-device distances;
- leaf/jaw positions каждого control point;
- SAD.

Validator для Versa HD/Agility ожидает 80 пар и nominal projected leaf width
5 mm с общим span 400 mm. Это delivery-coordinate check, а не full-head model.

## Energy mapping

Наличие одновременно 6 MV и 6 FFF (аналогично 10 MV / 10 FFF) означает, что
`NominalBeamEnergy=6` само по себе недостаточно.

Mapping использует:

1. NominalBeamEnergy;
2. PrimaryFluenceModeSequence / FluenceMode;
3. FluenceModeID;
4. при неоднозначности — только явный override.

Молчаливо угадывать FFF по имени beam запрещено.

## Reference Monte Carlo backend

Первым физическим эталоном остаётся EGSnrc/BEAMnrc + DOSXYZnrc.

Открытая литература показывает практически тот же workflow, который нужен
проекту:

- full Versa HD head model в BEAMnrc;
- SYNCMLCE / SYNCJAWS для dynamic delivery;
- DICOM RT Plan -> sequence files;
- patient CT -> DOSXYZnrc;
- сравнение с Monaco dose-to-medium.

Ключевая работа:
Paschal et al., 2022, *Monte Carlo modeling of the Elekta Versa HD and patient
dose calculation with EGSnrc/BEAMnrc*, JACMP, DOI 10.1002/acm2.13715.

Для Agility опубликованы модели с 80 leaf pairs, 5-mm projected leaves и
валидацией transmission / tongue-and-groove / IMRT / VMAT. Однако точные
manufacturer head dimensions в части работ были proprietary, поэтому
литературные числа нельзя без commissioning объявлять геометрией нашей машины.

## Практический вывод

M3 делится на два слоя:

1. **Delivery model** — DICOM geometry, Agility motion, jaws, gantry.
2. **Physical head model** — target/source, primary collimation, FF/FFF,
   monitor/head structures, physical Agility leaf body/materials.

Первый слой уже можно строить полностью. Второй должен быть получен из
документации/измерений и независимо валидирован PDD, profiles, output factors
и MLC-specific tests.
