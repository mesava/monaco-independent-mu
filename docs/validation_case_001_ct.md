# Validation case 001 — CT / RTSTRUCT

> **Статус: НЕ ВЕРИФИЦИРОВАНО.**
>
> Этот файл ранее содержал численные характеристики пациентской CT/RTSTRUCT,
> но в текущей воспроизводимой цепочке нет доступного raw DICOM набора, на
> котором эти числа можно заново получить и проверить. Поэтому прежние
> численные значения не считаются результатом проекта и не должны
> использоваться для физического или клинического вывода.

## Что уже реализовано в коде

Для реального validation case проект теперь умеет воспроизводимо:

- найти DICOM case и пройти цепочку RTPLAN → RTSTRUCT → referenced CT series;
- загрузить выбранный CT SeriesInstanceUID, не смешивая соседние CT series;
- сортировать CT по ImageOrientationPatient / ImagePositionPatient;
- применять RescaleSlope / RescaleIntercept;
- проверять FrameOfReferenceUID, геометрию и равномерность slice spacing;
- перечислять ROI без сохранения PatientName/PatientID;
- растеризовать CLOSED_PLANAR / CLOSEDPLANAR_XOR RTSTRUCT ROI;
- применить scanner-specific DRT120kV HU → RED calibration;
- построить RED → mass density → material model;
- явно учитывать HU вне measured calibration range;
- сформировать de-identified .egsphant, PEGSless media и JSON summary.

## Что требуется для закрытия case 001

Нужен повторный запуск на исходном пациентском DICOM наборе, доступном
вычислительной среде. После этого в этот файл можно записывать только
результаты, автоматически полученные текущей версией кода:

1. CT series UID / FrameOfReference consistency;
2. shape, spacing, orientation, KVP;
3. ROI name, contour count, rasterized volume;
4. HU min/max и out-of-calibration counts внутри Patient ROI;
5. mass-density range и material statistics;
6. hashes сгенерированных patient-model artifacts;
7. git commit, которым выполнен расчёт.

До такого повторного запуска **M1 остаётся реализованным на уровне кода и
synthetic integration tests, но не закрытым как real-patient validation**.

## Политика данных пациента

Исходные DICOM и персональные идентификаторы не коммитятся в репозиторий.
В Git допускаются только обезличенные конфигурации и агрегированные результаты,
которые можно воспроизвести из локального validation dataset.
