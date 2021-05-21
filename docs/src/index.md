# Replication of "Public Debt, Interest Rates, and Negative Shocks" (Evans, R. 2020) 

> This replication study was part of our evaluation for the course [Numerical Methods](https://floswald.github.io/NumericalMethods/) at SciencesPo Paris in Spring 2021
> 
> The functions used to replicate this paper are:

```@autodocs
Modules = [Evans2020]
```

## Case number 1 : &epsilon; = 1 and &mu; is constant

```julia
julia> using Evans2020

julia> 
```

end


# Our Replication of The Cash Paradox (Jiang & Shao, 2019)

> This replication study was part of our evaluation for the course [Numerical Methods](https://floswald.github.io/NumericalMethods/) at SciencesPo Paris in Spring 2021

The functions used to replicate this paper are:

```@autodocs
Modules = [CashParadox]
```

## Replication of model predictions: Figures 5a-5d

### Figure 5a: AUSTRALIA

In order to create this figure, type

```julia
julia> using Cash Paradox

julia> Fig5(2)
```

![fig5aus](./assets/fig5aus.png)

### Figure 5b: CANADA

In order to create this figure, type

```julia
julia> using Cash Paradox

julia> Fig5(0)
```

![fig5can](./assets/fig5can.png)

### Figure 5c:  UK

In order to create this figure, type

```julia
julia> using Cash Paradox

julia> Fig5(3)
```

![fig5uk](./assets/fig5uk.png)

### Figure 5d: US

In order to create this figure, type

```julia
julia> using Cash Paradox

julia> Fig5(1)
```
![fig5us](./assets/fig5us.png)


## Replication of model predictions: Figures A2a-A2d


### Figure A2a: AUSTRALIA

In order to create this figure, type

```julia
julia> using Cash Paradox

julia> FigA2(2)
```

![figA2aus](./assets/figA2aus.png)

### Figure A2b: CANADA

In order to create this figure, type

```julia
julia> using Cash Paradox

julia> FigA2(0)
```

![figA2can](./assets/figA2can.png)

### Figure A2c:  UK

In order to create this figure, type

```julia
julia> using Cash Paradox

julia> FigA2(3)
```
![figA2uk](./assets/figA2uk.png)

### Figure A2d: US

In order to create this figure, type

```julia
julia> using Cash Paradox

julia> FigA2(1)
```

![figA2us](./assets/figA2us.png)



## Replication of regime changes over time: Figure A3

In order to create this figure, type

```julia
julia> using Cash Paradox

julia> FigA3()
```

![figA3](./assets/figA3.png)

## Replication of The Value of ATM withdrawals over CIC: Figure A4

In order to create this figure, type

```julia
julia> using Cash Paradox

julia> FigA4()
```

![figA4](./assets/figA4.png)


## Replication of Cash receipts from circulation in the Federal Reserve Banks: Figure A5

In order to create this figure, type

```julia
julia> using Cash Paradox

julia> FigA5()
```

![figA5](./assets/figA5.png)


## Replication of Different measures of nominal interest rates: Figure D1

In order to create this figure, type

```julia
julia> using Cash Paradox

julia> FigD1()
```

![figD1](./assets/figD1.png)

## Replication of CIC/GDP with different interest rate specifications: Figure D2

In order to create this figure, type

```julia
julia> using Cash Paradox

julia> FigD2()
```
![figD2](./assets/figD2.png)

## Replication of Table 1 and Table 2: Calibration results and Cash shares relative to credit 

Most of results in this table are very close to those from Jiang & Shao (2019). However, we encounter problems replicating the result for the NCF model with UK data.

In order to create a dataframe containing this results type

```julia
julia> using Cash Paradox

julia> Table1()
```

![figG1](./assets/t1.png)

## Replication of Table D.1: Parameter values and Table D.2: Model performance comparison


In order to create a dataframe containing this results type

```julia
julia> using Cash Paradox

julia> Table2()
```

![figD1D2](./assets/t2.png)


