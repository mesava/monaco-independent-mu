# Validation case 001 — проверенный DICOM preflight

> Источник: полный обезличиваемый DICOM export, повторно предоставленный для
> проекта 2026-10-07.
>
> В репозиторий не коммитятся исходные DICOM, PatientName, PatientID или
> исходные имена файлов. Ниже записаны только физико-технические характеристики,
> которые были непосредственно извлечены из предоставленного архива.

## Комплект

Архив содержит **190 DICOM objects**:

- CT: **187**;
- RTSTRUCT: **1**;
- RTPLAN: **1**;
- RTDOSE: **1**.

## CT

Проверено:

- одна CT SeriesInstanceUID;
- один FrameOfReferenceUID;
- **187 slices**;
- matrix: **512 × 512**;
- PixelSpacing: **1.5625 × 1.5625 mm**;
- spacing along slice normal: **2.5 mm**, uniform;
- orientation: axial;
- PatientPosition: **HFS**;
- KVP: **120 kV**;
- RescaleSlope: **1**;
- RescaleIntercept: **−1024 HU**;
- coordinate along slice normal: **−212.5 … +252.5 mm**.

Таким образом этот case совместим с текущей scanner-specific
DICOM3.DRT120kV calibration по KVP и геометрически подходит для HFS
Source21 transform.

## RTSTRUCT

RTSTRUCT:

- использует тот же FrameOfReferenceUID, что CT;
- ссылается на CT Series из архива;
- reference list содержит все **187 CT SOP instances**.

В StructureSetROISequence находятся 14 ROI:

- Breast_L;
- Esophagus;
- Glnd_Thyroid;
- Heart;
- Humer_R;
- Liver;
- Lung_L;
- Lung_R;
- Lungs;
- SpinalCord;
- Trachea;
- Breast_R;
- CTV_LN;
- PTV_50Gy.

Все проанализированные contours имеют geometric type CLOSED_PLANAR.

### Важный blocker M1

В этом RTSTRUCT **нет отдельного external / BODY / patient ROI**.

Текущий production-oriented patient-model workflow требует явную внешнюю
маску, потому что voxels вне неё переводятся в DryAir, а couch/PPS должен
оставаться отдельным геометрическим слоем.

Поэтому нельзя молча:

- считать весь CT FOV пациентом;
- использовать PTV/грудь вместо внешнего контура;
- автоматически включать CT table/couch в Patient material LUT.

Для этого case M1 должен либо получить independently-derived external contour,
либо отдельный проверенный full-FOV/PPS strategy.

## RTPLAN

План содержит:

- planned fractions: **25**;
- treatment beams: **8**;
- RadiationType: **PHOTON**;
- BeamType: **DYNAMIC**;
- nominal energy: **6 MV**;
- Primary Fluence Mode: **STANDARD**;
- SAD: **1000 mm**;
- общий isocenter первого control point:
  **(−71.3, −84.5, +100.0) mm** в DICOM patient coordinates.

Это **не VMAT**: gantry angle внутри каждого beam постоянен, а движение
происходит в MLC. Для данного validation case delivery class — **dynamic MLC
IMRT**.

### Agility delivery geometry

Для каждого beam:

- beam-limiting MLC device: MLCX;
- leaf pairs: **80**;
- LeafPositionBoundaries: **−200 … +200 mm**;
- все projected leaf widths: **5 mm**;
- SourceToBeamLimitingDeviceDistance: **349 mm**.

Для Y jaws:

- device: ASYMY;
- one jaw pair;
- SourceToBeamLimitingDeviceDistance: **470 mm**.

Это реальным DICOM export подтверждает ожидаемую delivery geometry Agility:
80 × 5-mm projected leaves over a 400-mm span.

### Beam MU и control points

| Beam | MU | CP | Gantry | Collimator |
|---:|---:|---:|---:|---:|
| 1 | 127.458344 | 30 | 51° | 348° |
| 2 | 70.082390 | 30 | 226° | 12° |
| 3 | 58.973721 | 29 | 206° | 10° |
| 8 | 74.168732 | 29 | 240° | 10° |
| 4 | 45.371567 | 18 | 75° | 348° |
| 5 | 127.043335 | 30 | 0° | 0° |
| 6 | 98.906357 | 32 | 270° | 0° |
| 7 | 92.862953 | 33 | 245° | 10° |

Суммарно: **694.867399 MU per fraction**.

Couch для beams находится при 0°. CumulativeMetersetWeight идёт от 0 до 1,
а MLC states меняются между control points.

## RTDOSE

В архиве присутствует один RTDOSE:

- DoseUnits: **GY**;
- DoseType: **PHYSICAL**;
- DoseSummationType: **PLAN**;
- FrameOfReferenceUID совпадает с CT/RTSTRUCT/RTPLAN.

То есть этот export пригоден для **plan-level M6 3D comparison**.

Но он **не содержит per-beam RTDOSE objects**, поэтому сам по себе не
позволяет выполнить patient-specific M5 equivalent-MU numerator для каждого
beam. Для M5 позднее понадобится экспорт TPS dose per beam или эквивалентный
надёжный источник per-beam TPS dose.

## Что этот case уже закрыл

Реальный Monaco export подтвердил:

1. DICOM case discovery должен выбрать одну CT series из RTSTRUCT reference;
2. HFS / 120-kV transport preflight подходит;
3. MLC device действительно экспортируется как MLCX, 80 pairs;
4. Agility boundaries действительно имеют 5-mm spacing;
5. план может быть dynamic MLC IMRT без gantry rotation;
6. один plan RTDOSE достаточен для M6, но не для M5 per-beam comparison.

## Следующий шаг

1. реализовать независимый external-contour strategy как отдельный,
   валидируемый слой, поскольку patient/BODY отсутствует;
2. после появления patient mask прогнать реальный CT:
   HU → RED → rho → material → .egsphant;
3. сформировать machine transport input для этого 8-field DMLC IMRT;
4. сравнить independent MC с plan RTDOSE;
5. отдельно получить per-beam TPS dose для M5.


## Предварительный CT-derived external-mask прогон

Исходный архив повторно распакован локально и current research algorithm
воспроизведён непосредственно на pixel data CT без коммита DICOM в Git.

Конфигурация:

- seed: общий treatment isocenter **(-71.3, -84.5, +100.0) mm**;
- threshold: **-500 HU**;
- 3-D closing iterations: **1**;
- 6-connected component;
- enclosed holes filled per axial slice.

Получено:

- CT observed HU range: **-3024 ... +3071 HU**;
- HU в voxel ближайшем к isocenter: **-821 HU**;
- поэтому seed находится в low-density lung region и не принадлежит
  thresholded component;
- nearest-tissue snap distance: **2.210 mm**;
- selected component before per-slice hole fill: **5,898,483 voxels**;
- after fill: **6,435,424 voxels**;
- volume: **39,278.71 cm³**;
- mask does **not** touch CT volume border;
- mask z extent: **-210 ... +250 mm**;
- HU range inside the filled mask: **-1024 ... +3071 HU**.

Последний пункт важен: scanner calibration DRT120kV имеет measured range
**-1000 ... +2009 HU**, поэтому real patient mask содержит voxels за обоими
концами calibration table. Production code не должен молча это скрывать:
counts ниже/выше calibration range выводятся отдельно, а endpoint clipping
разрешается только явной policy.

### Интерпретация

Seed-snap механизм для этого case оказался действительно необходимым:
treatment isocenter расположен в лёгком и имеет HU ниже -500.

Полученный объём около 39.3 L физически правдоподобен для длинного breast/chest
CT с руками, но **сам по себе не доказывает отсутствие CT couch/support в
маске**. Поэтому результат пока имеет статус:

    RESEARCH_MASK_DERIVED_NOT_VALIDATED

До записи patient .egsphant как validation artifact необходимо:

1. прогнать threshold sensitivity;
2. проверить posterior support contamination;
3. зафиксировать counts HU < -1000 и HU > 2009 внутри mask;
4. только затем выбрать явную out-of-range policy и material-mixture
   discretisation.
