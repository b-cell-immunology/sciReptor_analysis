# Name:			V-segment_association_heatmap.R
# Author:		Christian Busse
# Maintainer:	Christian Busse
# Version:		0.1.4
# Last change:	2016-02-16
# License:		AGPLv3
# Requires:		Database with a "derived_segments_view" table
# Description:	This script produces heatmaps of associated heavy/light V-segments.
#
library(RMySQL)

# Config parameters
#
config.mysql.group <- "mysql_igdb"
config.output.path <- Sys.getenv("HOME")
config.database <- "omega_w05"
config.select <- "tissue='spleen' AND (population='MN' OR population= 'B')"
#config.select <- "tissue='spleen' AND population='MN'"


# If available, load KWallet authentication library
#
module.auth.kwallet.path <- file.path(
	Sys.getenv("HOME"),
	"development/lab_public/lib_authentication_kwallet.R"
)
if (file.exists(module.auth.kwallet.path)) {
	source(module.auth.kwallet.path)
}


# Function to query a unified list of associated V segments from a database. The function will automatically select
# the existing and functional light chain locus for the 'light_*' columns which are returned. Note that in case of
# both kappa and lambda being present and functional, kappa will be given precedence. Further note that the optional
# 'select.current' parameter uses the column names of the 'derived_segments_view' tables not the ones returned by the
# function.
#
func.get.list.unified <- function(connection.current, db.current, select.current) {
	if(missing(select.current)) {
		select.current <- "TRUE"
	}
	temp.df <- dbGetQuery(
		connection.current,
		paste(
			"SELECT ",
				"CAST(event_id AS UNSIGNED INTEGER) AS event_id, ",
				"donor_identifier, ",
				"tissue, ",
				"population, ",
				"igh_segment_v AS heavy_segment_v, ",
				"igh_segment_d AS heavy_segment_d, ",
				"igh_segment_j AS heavy_segment_j, ",
				"igh_cdr3 AS heavy_cdr3, ",
				"igh_segment_c as heavy_constant, ",
				"igk_segment_v AS light_segment_v, ",
				"igk_segment_j AS light_segment_j, ",
				"igk_cdr3 AS light_cdr3, ",
				"igk_segment_c as light_constant, ",
				"CONCAT(igh_segment_v,'_',igk_segment_v) AS vv_associated ",
			"FROM ", db.current, ".derived_segment_association ",
			"WHERE igk_segment_v IS NOT NULL AND (", select.current,") ",
			"UNION ",
			"SELECT ",
				"CAST(event_id AS UNSIGNED INTEGER) AS event_id, ",
				"donor_identifier, ",
				"tissue, ",
				"population, ",
				"igh_segment_v AS heavy_segment_v, ",
				"igh_segment_d AS heavy_segment_d, ",
				"igh_segment_j AS heavy_segment_j, ",
				"igh_cdr3 AS heavy_cdr3, ",
				"igh_segment_c as heavy_constant, ",
				"igl_segment_v AS light_segment_v, ",
				"igl_segment_j AS light_segment_j, ",
				"igl_cdr3 AS light_cdr3, ",
				"igl_segment_c as light_constant, ",
				"CONCAT(igh_segment_v,'_',igl_segment_v) AS vv_associated ",
			"FROM ", db.current,".derived_segment_association ",
			"WHERE igk_segment_v IS NULL AND igl_segment_v IS NOT NULL AND (", select.current,");",
			sep=""
		)
	)
	return(temp.df)
}

# Function to convert a vector of associated elements (one-dimensional, elements are separated by single character (default '_'))
# into a two-dimensional matrix, counting the number of observed associations.
#
func.gen.assoc.matrix <- function(vector.string.A, vector.string.B) {
	if (length(vector.string.A) != length(vector.string.B)) {
		stop("Input strings have non-equal length")
	}
	if (length(c(grep("\t", vector.string.A), grep("\t", vector.string.B)))) {
		stop("Input strings contain tabs which are used as internal separators")
	}

	vector.string.uniq.A <- unique(vector.string.A)
	vector.string.uniq.B <- unique(vector.string.B)
	matrix.output <- matrix(
		rep(0, length(vector.string.uniq.A) * length(vector.string.uniq.B)),
		nrow=length(vector.string.uniq.A),
		dimnames=list(vector.string.uniq.A, vector.string.uniq.B)
	)

	vector.strings.freq <- table(paste(vector.string.A,vector.string.B,sep="\t"))
	matrix.assoc.counts <- cbind(
		matrix(unlist(strsplit(names(vector.strings.freq),"\t")),ncol=2,byrow=TRUE),
		unname(vector.strings.freq)
	)
	matrix.output[matrix.assoc.counts[,c(1,2)]]<-as.integer(matrix.assoc.counts[,3])
	return(matrix.output)
}


# Connect to database. If authentication via KWallet is available, it will be used, otherwise
# the password has to be provided as "password" key in the .my.cnf config file.
#
if (exists("auth.kwallet.avail") && auth.kwallet.avail) {
	connection.mysql <- dbConnect(
		MySQL(),
		group = config.mysql.group,
		password = getKWalletPassword("kdewallet", "MySQL", generateMySQLKWalletKey(config.mysql.group))
		)
} else {
	connection.mysql <- dbConnect(
		MySQL(),
		group = config.mysql.group,
	)
}


# == Query data and create various association matrices ==
#
df.associated <- func.get.list.unified(connection.mysql, config.database, config.select)
observations.total <- nrow(df.associated)

# calculate the expected frequency based on random association, as the outer product of the heavy and light "count"
# vectors. Since the sum of the individual vectors is n (total number of events) the sum ob the resulting "count" matrix
# is n^2. This is then normalized in the "freq" matrix. The "atom" is the fraction 1/n^2 which in the "freq" matrix is then
# number that corresponds to 1 in the count matrix. As such, all numbers in the matrix are an integer multiple of it, however
# it is not necessarily (but likely) the largest common divisor (as far as such a thing exists for floats). It will be used
# later on in operations to avoid divisions by zero, without altering the underlying data too much or at least in a predictable
# manner.
#
matrix.count.expected <- outer(
	rev(sort(table(factor(df.associated[,"heavy_segment_v"])))),
	rev(sort(table(factor(df.associated[,"light_segment_v"]))))
)
matrix.freq.expected <- matrix.count.expected / sum(matrix.count.expected)
matrix.atom <- 1 / sum(matrix.count.expected)

# perform absolute and relative counting of VV segment associations
#
matrix.count.observed <- func.gen.assoc.matrix(df.associated[,"heavy_segment_v"],df.associated[,"light_segment_v"])
matrix.count.observed <- matrix.count.observed[rownames(matrix.count.expected), colnames(matrix.count.expected)]
matrix.freq.observed <- matrix.count.observed / sum(matrix.count.observed)

# Alternative ratio calculation, add atomic value to complete matrix before division. This is a comparable strategy to log(N+1)
# in terms of avoiding log(0), but has a smaller effect on associations with an higher expected frequency (since it adds before the division).
#
matrix.freq.ratio <- (matrix.freq.observed+matrix.atom) / matrix.freq.expected
matrix.freq.ratio.log <- log(matrix.freq.ratio,2)

# Calculate likelihood that NO events of a given combination would be sampled, given the actual sampling depth (n) and the expected
# likelihood based on the one-dimensional (heavy or light only) frequencies.
#
matrix.likelihood.miss <- (1-matrix.freq.expected)^(observations.total)

# Filter out (i.e. set to 0) all events where that would not be expected to be found at a sampling depth 3-fold higher than actual
#
matrix.freq.ratio.log[which(matrix.likelihood.miss > 0.95) ]<-0


# == Define internal output parameters, open and format output. ==
#
plot.margin.bottom<-0.25
plot.margin.left.map<-0.875
plot.margin.left.legend<-0.125
plot.margin.right<-0.125
plot.margin.top<-1.00

pdf(file=file.path(config.output.path, "vv_associated_heatmap_1_uf.pdf"), paper="A4r", width=11.7, height=8.27)

layout(matrix(c(1,2),nrow=1),width=c(15,1))

# == Output absolute counts heatmap ==
#

# Transpose observed counts matrix and invert row order for output
#
matrix.output <- t(matrix.count.observed)[,nrow(matrix.count.observed):1]

heatmap.breaks<-c(-1:max(matrix.output))
palette.current <-rev(c(
	gray((0:max(matrix.output))/max(matrix.output))
))

par(mai=c(plot.margin.bottom,plot.margin.left.map,plot.margin.top,plot.margin.right),ps=7,omi=c(0,0,0,0.25),las=2)
image(matrix.output, col=palette.current, breaks=heatmap.breaks, axes=FALSE)
box()
axis(2,seq(from=0,to=1,length.out=ncol(matrix.output)),colnames(matrix.output))
axis(3,seq(from=0,to=1,length.out=nrow(matrix.output)),rownames(matrix.output))
title(
	main=paste(config.database, " (n=",nrow(df.associated),")",sep=""),
	outer=FALSE,
	font.main=2,
	cex.main=2,
	line=4
)
title(
	sub=paste("selected for \"",config.select,"\"",sep=""),
	line=0,
	outer=FALSE,
	cex.sub=1.25
)

par(mai=c(plot.margin.bottom,plot.margin.left.legend,plot.margin.top,plot.margin.right),ps=9,font=2)
image(matrix(heatmap.breaks[-1],nrow=1),col=palette.current,breaks=heatmap.breaks,axes=FALSE)
box()
title(main="count",outer=FALSE,font.main=2,line=1)
axis(4,(0:max(matrix.output))/max(matrix.output),0:max(matrix.output))


# == Output expected likelihood heatmap ==
#

# Transpose log ratio matrix and invert row order for output
#
matrix.output <- t(matrix.freq.ratio.log)[,nrow(matrix.freq.ratio.log):1]

# heatmap.breaks<-rev(c(308,log(2^(5:-5)+1,2),0)) # The max of 308 is based on 64-bit float maxima
# palette.current <-rev(c(
# 	rgb(1,seq(from=0,to=0.81,length.out=5),seq(from=0,to=0.81,length.out=5)),
# 	gray(0.75),
# 	rev(rgb(seq(from=0,to=0.81,length.out=5),seq(from=0,to=0.81,length.out=5),1)),
# 	gray(1)
# ))

heatmap.breaks<-c(-12,-10,-8,-6,-4,-2,-1,1,2,4,6,8,10,12)
palette.current <-rev(c(
	rgb(1,seq(from=0,to=0.81,length.out=6),seq(from=0,to=0.81,length.out=6)),
	gray(1),
	rev(rgb(seq(from=0,to=0.81,length.out=6),seq(from=0,to=0.81,length.out=6),1))
))

par(mai=c(plot.margin.bottom,plot.margin.left.map,plot.margin.top,plot.margin.right),ps=7,omi=c(0,0,0,0.25),las=2)
image(matrix.output, col=palette.current, breaks=heatmap.breaks, axes=FALSE)
box()
axis(2,seq(from=0,to=1,length.out=ncol(matrix.output)),colnames(matrix.output))
axis(3,seq(from=0,to=1,length.out=nrow(matrix.output)),rownames(matrix.output))
title(
	main=paste(config.database, " (n=",nrow(df.associated),")",sep=""),
	outer=FALSE,
	font.main=2,
	cex.main=2,
	line=4
)
title(
	sub=paste("selected for \"",config.select,"\"",sep=""),
	line=0,
	outer=FALSE,
	cex.sub=1.25
)

par(mai=c(plot.margin.bottom,plot.margin.left.legend,plot.margin.top,plot.margin.right),ps=9,font=2)
image(matrix(heatmap.breaks[2:13],nrow=1),col=palette.current,breaks=heatmap.breaks,axes=FALSE)
box()
title(main=" log2(x+1)",outer=FALSE,font.main=2,line=1)
axis(4,seq(from=0,to=1,length.out=12),rev(c(as.character(5:-5),"n/a")))


# == Clean up ==
#

if(any(names(dev.list())=="pdf")){
	dev.off(unname(dev.list()[names(dev.list())=="pdf"]))
}

dbDisconnect(connection.mysql)
