# DOSXYZnrc 3ddose

Добавлен собственный parser/writer стандартного dense ASCII `.3ddose`.

Формат DOSXYZnrc:

1. `nx ny nz`;
2. `nx+1` x-boundaries, cm;
3. `ny+1` y-boundaries, cm;
4. `nz+1` z-boundaries, cm;
5. `nx*ny*nz` dose values;
6. столько же relative statistical uncertainties.

В файле x меняется быстрее y, а y быстрее z. Внутри проекта массив хранится
как:

```text
dose[z, y, x]
```

## Важная политика единиц

Parser не называет сырое значение `Gy/MU`.

Выход DOSXYZnrc сначала является результатом конкретной MC normalization.
Клиническая абсолютная доза появляется только после применения независимо
полученной calibration factor.

Поэтому `ThreeDDose.scaled()` выполняет явное масштабирование dose, не меняя
relative statistical uncertainty.

Это позволит в дальнейшем отдельно хранить:

- raw MC dose/history;
- independent absolute calibration factor;
- final Gy/MU;
- beam MU;
- final beam/plan Gy.

## Источник формата

NRC DOSXYZnrc User Manual, section 12.1, Format of .3ddose.
