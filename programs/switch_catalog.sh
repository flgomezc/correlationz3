#

clear
if [ -d ~/correlationz3/data_small ];   then
    echo "#######   Switch to TEST data"
    mv ~/correlationz3/data       ~/correlationz3/data_full
    mv ~/correlationz3/data_small ~/correlationz3/data
else
    echo "#######   Switch to FULL data"
    mv ~/correlationz3/data      ~/correlationz3/data_small
    mv ~/correlationz3/data_full ~/correlationz3/data
fi
echo
