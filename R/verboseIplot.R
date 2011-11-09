verboseIplot<-function (
x, 
y, 
xlim=NA,
ylim=NA,
nBinsX=150, 
nBinsY=150,
ztransf = function(x){x}, 
gamma=1,
sample = NULL, 
corFnc = "cor", 
corOptions = "use = 'p'", 
main = "", 
xlab = NA, 
ylab = NA, 
cex = 1, 
cex.axis = 1.5, 
cex.lab = 1.5, 
cex.main = 1.5, 
abline = FALSE, 
abline.color = 1, 
abline.lty = 1, 
corLabel = corFnc,
...) 
{	
	if (is.na(xlab))         
		xlab = deparse(substitute(x))    
	if (is.na(ylab))         
		ylab = deparse(substitute(y))           

	x = as.numeric(as.character(x))    
	y = as.numeric(as.character(y))    
	xy <- data.frame(x,y)    
	xy=xy[!is.na(x)&!is.na(y),]        

	if (sum(is.na(xlim))!=0)         
		xlim=c(min(xy[,1])-10^-10*diff(range(xy[,1])),max(xy[,1]))    
	if (sum(is.na(ylim))!=0)         
		ylim=c(min(xy[,2])-10^-10*diff(range(xy[,2])),max(xy[,2]))                    

	corExpr = parse(text = paste(corFnc, "(x, y ", prepComma(corOptions),")"))    
	cor = signif(eval(corExpr), 2)    
	corp = signif(corPvalueStudent(cor, sum(is.finite(x) & is.finite(y))),2)    
	if (corp < 10^(-200))         
		corp = "<1e-200"    
	else corp = paste("=", corp, sep = "")
   
	 resid=lm(y~x)$residuals
	MSE=round(mean(resid^2),2)
	if (!is.na(corLabel)) {        
		mainX = paste(main, " ", corLabel, "=", cor, " MSE = ", MSE,sep = "")        }    
	else mainX = main        

	if (!is.null(sample)) {        
		if (length(sample) == 1) {            
			sample = sample(length(x), sample)        }        
	xy=xy[sample,]           }                

	sx <- seq(xlim[1], xlim[2], by = diff(xlim)/nBinsX)    
	sy <- seq(ylim[1], ylim[2], by = diff(ylim)/nBinsY)         
	den <- ztransf(table(cut(xy[, 1], breaks =  sx), cut(xy[, 2], breaks = sy)))        

	lsx <- length(sx)    
	lsy <- length(sy)    
	xx <- 0.5 * (sx[-1] + sx[-lsx])    
	yy <- 0.5 * (sy[-1] + sy[-lsy])      
	
	whiteBlueGreenRedBlack=function (n){      
		quarter = as.integer(n/5)      
		red=  c(seq(from=1,to=0,length.out=quarter)^(1/gamma),seq(from=0,to=0,length.out=quarter)^(1/gamma),seq(from=0,to=1,length.out=quarter)^(1/gamma),seq(from=1,to=1,length.out=quarter)^(1/gamma),seq(from=1,to=0,length.out=quarter)^(1/gamma))       		
		green=c(seq(from=1,to=1,length.out=quarter)^(1/gamma),seq(from=1,to=1,length.out=quarter)^(1/gamma),seq(from=1,to=1,length.out=quarter)^(1/gamma),seq(from=1,to=0,length.out=quarter)^(1/gamma),seq(from=0,to=0,length.out=quarter)^(1/gamma))      		
		blue= c(seq(from=1,to=1,length.out=quarter)^(1/gamma),seq(from=1,to=0,length.out=quarter)^(1/gamma),seq(from=0,to=0,length.out=quarter)^(1/gamma),seq(from=0,to=0,length.out=quarter)^(1/gamma),seq(from=0,to=0,length.out=quarter)^(1/gamma))      
		col = rgb(red, green, blue, maxColorValue = 1)      
		col    }        


	image(x = xx, y = yy, den,  xaxs = "r", yaxs = "r", xlab = xlab, ylab = ylab, cex = cex,     main=mainX, cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main, col=whiteBlueGreenRedBlack(50))           

	if (abline) {        	
	fit = lm(y ~ x)        
	abline(reg = fit, col = abline.color, lty = abline.lty)    }    
	invisible(sample)
}
