# External patient mask strategy

## Почему слой вынесен отдельно

Для patient Monte Carlo недостаточно просто загрузить всю CT matrix.

Вне пациента CT может содержать:

- air background;
- scanner padding;
- table/support structures;
- immobilization objects.

При этом Patient material LUT и PPS/support material model имеют разные
физические роли.

Поэтому external patient mask является отдельным входом patient model.

## Preferred production source

Предпочтительный путь:

1. RTSTRUCT содержит отдельный external/BODY/patient ROI;
2. ROI ссылается на тот же FrameOfReferenceUID и CT series;
3. контуры растеризуются на native CT grid;
4. voxels вне external ROI переводятся в DryAir;
5. couch/PPS моделируется отдельно.

## Validation case 001

В полученном реальном Monaco export отдельный external/BODY/patient ROI
отсутствует.

Нельзя молча заменить его:

- всем CT FOV;
- PTV;
- breast contour;
- объединением внутренних органов.

Поэтому добавлен только явно выбираемый research path.

## CT-derived research mask

Функция derive_external_mask_from_ct использует:

1. явный HU threshold;
2. optional 3-D binary closing;
3. connected-component labeling;
4. treatment isocenter как patient-space seed;
5. если isocenter находится в low-density cavity, разрешается привязка к
   ближайшему thresholded component только в пределах явного maximum distance;
6. enclosed holes заполняются по CT slices.

Это позволяет не выбирать largest component скрытой эвристикой и может
отбросить отдельный CT support object.

## Ограничение

Если тело пациента физически соединяется с CT table/support на выбранном HU
threshold, connected-component segmentation может включить support вместе с
пациентом.

Поэтому CT-derived mask имеет статус:

    RESEARCH_ONLY

и не является clinical default.

## Что нужно для его валидации

Нужен набор CT, где одновременно доступны:

- trusted external contour;
- original CT;
- treatment isocenter.

Для серии пациентов следует сравнить CT-derived и reference masks:

- Dice;
- surface Dice;
- Hausdorff / 95th percentile surface distance;
- volume difference;
- posterior support contamination;
- differences в HU/RED/material map;
- dose sensitivity в independent MC.

Отдельно нужны сложные случаи:

- breast board;
- vacuum bag;
- arms close to body;
- lung / large internal air cavities;
- truncation at CT FOV;
- metal implants;
- very low-density anatomy.

До такой validation series реальный case 001 может использовать этот path
только как исследовательский тест patient-model pipeline.
