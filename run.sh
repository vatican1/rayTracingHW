#!/bin/bash

# Проверка, что передано два аргумента
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <arg1> <arg2>"
    exit 1
fi

# Сохраняем аргументы
ARG1=$1
ARG2=$2

# Запускаем вашу собранную программу с переданными аргументами
./build/RayTracer "$ARG1" "$ARG2"
