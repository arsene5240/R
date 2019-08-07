## match.call() 用法：參數匹配
f <- function(x, y, ...) match.call()
f(y = 1, 2, z = 3, 4)
f(x = 2, y = 1, z = 3, 4)

## environment() 用法：環境變數
f <- function(f_x){
  g <- function(g_x){
    print("Inside g")
    print(environment())
    print(ls())
    message("......")
    }
  g(5)
  print("Inside f")
  print(environment())
  print(ls())
}
f(3)

## match.call(expand.dots = FALSE) 用法：使得所有 ... 的參數標記為 ... 的單個參數
f <- function(x, y, ...) match.call(expand.dots = FALSE)
f(y = 1, 2, z = 3, 4)
f(x = 2, y = 1, ... = list(z = 3, 4))

## match() 用法：類似於 %in%
print(match(c(5, 2, 4), c(2, 7, 5, 3), nomatch = 123))

## 1L 與 1 的差別在變數類型不同
class(5L)
class(5)
x <- 1:100
typeof(x) # integer
y <- x+1
typeof(y) # double, twice the memory size
object.size(y) # 840 bytes (on win64) 
z <- x+1L
typeof(z) # still integer
object.size(z) # 440 bytes (on win64) 

## eval() 用法：可以將字串執行成 R 的指令 (evaluates an expression)
eval(1+1) # 1+1 is expression
eval("2+2") # "2+2" is string

## parse() 用法：可以將字串轉化為 expression (change the string into an expression)
eval(parse(text = 1+1)) # 1+1 is expression
parse(text = "2+2") # "2+2" is expression

##
