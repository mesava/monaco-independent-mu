# Архитектура

## Принцип

Проект разделяет четыре независимых слоя:

1. **DICOM / геометрия** — только чтение и нормализация данных.
2. **Patient model** — CT → HU → RED → density → material.
3. **Beam/transport model** — Versa HD / Agility / Monte Carlo.
4. **QA** — MU-check и 3D comparison.

Ни один слой QA не должен менять исходную MC-дозу или выполнять best-fit нормировку к Monaco.

## Основной поток данных

```text
RTPLAN ──────────────► Beam/CP model ─────┐
                                         │
CT ─► Patient model ─► MC voxel geometry ├─► Monte Carlo ─► 3D dose
                                         │
Machine commissioning ─► Beam model ─────┘

RTDOSE ──────────────────────────────────────────────► comparison
RTSTRUCT ───────────────────────────────────────────► ROI/DVH
```

## Независимость

Допустимо использовать из Monaco/DICOM:

- MU и delivery geometry;
- CT calibration конкретного сканера;
- измеренные commissioning data;
- исходный RT Dose как reference для сравнения.

Недопустимо использовать рассчитанную Monaco дозу как вход для независимого MC, кроме явно определённых сравнительных метрик.
