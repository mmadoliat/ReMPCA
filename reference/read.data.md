# Read ReMPCA Example Data

Loads one of the datasets included with the ReMPCA package.

## Usage

``` r
read.data(data = c("bike_day", "bike_hour", "energydata"))
```

## Arguments

- data:

  Character string specifying the dataset to load. Available options are
  `"bike_day"`, `"bike_hour"`, and `"energydata"`.

## Value

The requested dataset as a data frame.

## Examples

``` r
bike_day <- read.data("bike_day")
bike_hour <- read.data("bike_hour")
energy <- read.data("energydata")
```
