The scripts in wind computations:

- rotation_daily_v1: script to rotate wind speed components u,v, in along-shore and cross shore direction. In the script there is also computation of taux,tauy, tau_along and tau_cross (wind stress), for daily fields data.

- rotation_monthly_v1: to rotate vectors of monthly wind speed, and to compute the wind stress along shore and cross shore.

--------------------------------------------------------------

- ccmp_daily_products_script: modified from rotation_daily , is the script, which produce the new files with new variables, like the taux, tauy, ualong, ucross.

- ccmp_mothly_products_script: monthly version of the daily one.