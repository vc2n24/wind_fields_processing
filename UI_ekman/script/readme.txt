
- Only these two script are the ones to use for analysis

- UI_computation: script that calculate the UI from the masked_km data. It take the data from Z:\mpoc\rpdmoc\users\vc2n24\wind_fields\data_products\masked\masked_50km, and produce UI data in the folder: Z:\mpoc\rpdmoc\users\vc2n24\wind_fields\UI_ekman\UI_ekman_50km.

- UI_climatologies_monthly_method.m : script to calculate the climatologies from UI data. It calculate the daily and mnonthly climatologies, and anomalies. It calculate them by putting together all the measurment in the months over the years, finding the monthly climatology so 12 values, then padding the year and interpolating with smooth gaussian curve over the months.
it get the latitude band 29 - 25 climatoligies and anomalies, also the grid points ones.

-UI_climatologies_daily_method:  script to calculate the climatologies from UI data. It calculate the daily and monthly climatologies, and anomalies. It calculate them by putting together all the measurment in the months over the years, finding the daily climatology so 365. it get the latitude band 29 - 25 climatoligies and anomalies, also the grid points ones.


Th eother script have been used to run analysis and plot over old data, which were wrong due to wrong climatologies and rotation of winds. But some of them have nice plots and maps ideas.
