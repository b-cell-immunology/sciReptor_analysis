# Name:			V-segment_usage_bargraph.R
# Author:		Christian Busse
# Maintainer:	Christian Busse
# Version:		0.2.2
# Last change:	2015-08-11
# License:		AGPLv3
# Requires:		Database with a "derived_segments_view" table
# Description:	This script produces ordered bargraphs of V-segment usage of all populations and tissues in a
#				given data set.
#
library(RMySQL)

# Config parameters
#
config.mysql.group <- "mysql_igdb"
config.output.path <- Sys.getenv("HOME")
config.database <- "omega_w05"

# If available, load KWallet authentication library
#
module.auth.kwallet.path <- file.path(
	Sys.getenv("HOME"),
	"development/lab_public/lib_authentication_kwallet.R"
)
if (file.exists(module.auth.kwallet.path)) {
	source(module.auth.kwallet.path)
}


# The main plot function. Requires a two-column data frame with "segment" name and "counts", plus chain type and population strings (for output only)
#
make.plot <- function(usage.segment, chain.type, population) {
	par(las=2, oma=c(2,0,0.5,1))

	if(dim(usage.segment)[1] > 0){
		usage.max.abs <- max(usage.segment[,"counts"])
		usage.cumulate.rel <- cumsum(usage.segment[,"counts"]) / sum(usage.segment[,"counts"])
		index.half.use <- match(TRUE, usage.cumulate.rel >=0.5 )

		barplot(usage.segment[,"counts"],names.arg=usage.segment[,"segment"],space=0,axes=FALSE,ylab="count")
		lines(
			c(0, 1:dim(usage.segment)[1]-0.5, dim(usage.segment)[1]),
			c(0, usage.cumulate.rel * max(usage.segment[,"counts"]), max(usage.segment[,"counts"]))
		)
		segments(
			index.half.use-0.5, usage.segment[index.half.use,"counts"],
			index.half.use-0.5, usage.cumulate.rel[index.half.use]*usage.max.abs
		)
		segments(
			index.half.use-0.5, usage.cumulate.rel[index.half.use]*usage.max.abs,
			dim(usage.segment)[1]*1.1, usage.cumulate.rel[index.half.use]*usage.max.abs
		)
		segments(
			trunc((dim(usage.segment)[1]+1)/2,0)-0.5, usage.segment[trunc((dim(usage.segment)[1]+1)/2,0),"counts"],
			trunc((dim(usage.segment)[1]+1)/2,0)-0.5, usage.cumulate.rel[trunc((dim(usage.segment)[1]+1)/2,0)]*usage.max.abs
		)
		segments(
			trunc((dim(usage.segment)[1]+1)/2,0)-0.5, usage.cumulate.rel[trunc((dim(usage.segment)[1]+1)/2,0)]*usage.max.abs,
			dim(usage.segment)[1]*1.1, usage.cumulate.rel[trunc((dim(usage.segment)[1]+1)/2,0)]*usage.max.abs
		)

		axis(2)
		axis(4,seq(from=0,to=(max(usage.segment[,"counts"])),length.out=6),seq(from=0,to=100,length.out=6))
		mtext("cumulate [%]",side=4,las=0,line=2)
		mtext(
			paste(as.character(trunc(usage.cumulate.rel[trunc((dim(usage.segment)[1]+1)/2,0)]*100,1)),"%",sep=""),
			side=4, 
			at=usage.cumulate.rel[trunc((dim(usage.segment)[1]+1)/2,0)]*usage.max.abs,
			line=0.5,
			las=2
		)
		mtext(
			paste(as.character(trunc((index.half.use/dim(usage.segment)[1])*100,1)),"%",sep=""),
			side=1,
			at=index.half.use-0.5,
			line=-0.25,
			las=1
		)
	} else {
		plot.new()
		box()
		text(0.5, 0.5, labels=c("No data available"), adj=c(0.5,0.5))
	}
	title(
		paste(
			population,
			" ",
			chain.type,
			": Usage of ",
			dim(usage.segment)[1],
			" distinct segments (n=",
			sum(usage.segment[,"counts"]),
			")",
			sep=""
		)
	)
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


# Get all distinct populations and iterate through them
#
df_pop_tissue <- dbGetQuery(
	connection.mysql,
	paste(
		"SELECT DISTINCT population, tissue ",
		"FROM ", config.database, ".derived_segment_association;",
		sep=""
	)
)

list_pop_tissue<-lapply(c(1:dim(df_pop_tissue)[1]),function(x){df_pop_tissue[x,]})

sapply(list_pop_tissue, function(pop_tissue_current){ 
	dbGetQuery(
		connection.mysql,
		paste(
			"SET @POPULATION := '", pop_tissue_current[["population"]], "';",
			sep=""
		)
	)

	dbGetQuery(
		connection.mysql,
		paste(
			"SET @TISSUE := '", pop_tissue_current[["tissue"]], "';",
			sep=""
		)
	)

	temp_output_name <- paste(pop_tissue_current[["population"]], "from", pop_tissue_current[["tissue"]])
	temp_file_name <- paste(gsub("[[:blank:]]", "_", pop_tissue_current[["tissue"]]), "_", pop_tissue_current[["population"]], sep="")
	cat(paste("Processing population",temp_output_name, "\n"))

	usage.segment.heavy <- dbGetQuery(
		connection.mysql,
		paste(
			"SELECT ",
				"igh_segment_v AS segment, ",
				"count(event_id) AS counts ",
			"FROM ", config.database, ".derived_segment_association ",
			"WHERE population = @POPULATION ",
				"AND tissue = @TISSUE ",
			"GROUP BY segment ",
			"ORDER BY counts DESC;",
			sep=""
		)
	);

	usage.segment.kappa <-  dbGetQuery(
		connection.mysql,
		paste(
			"SELECT ",
				"igk_segment_v AS segment, ",
				"count(event_id) AS counts ",
			"FROM ", config.database, ".derived_segment_association ",
			"WHERE population = @POPULATION ",
				"AND tissue = @TISSUE ",
				"AND igk_segment_v IS NOT NULL ",
			"GROUP BY segment ",
			"ORDER BY counts DESC;",
			sep=""
		)
	);

	usage.segment.lambda <-  dbGetQuery(
		connection.mysql,
		paste(
			"SELECT ",
				"igl_segment_v AS segment, ",
				"count(event_id) AS counts ",
			"FROM ", config.database, ".derived_segment_association ",
			"WHERE population = @POPULATION ",
				"AND tissue = @TISSUE ",
				"AND igl_segment_v IS NOT NULL ",
			"GROUP BY segment ",
			"ORDER BY counts DESC;",
			sep=""
		)
	);

	pdf(file=paste("segment_usage_", config.database, "_", temp_file_name, ".pdf", sep=""), width=15, height=10, onefile=TRUE)
	make.plot(usage.segment.heavy, "heavy", paste(config.database, temp_output_name))
	make.plot(usage.segment.kappa, "kappa", paste(config.database, temp_output_name))
	make.plot(usage.segment.lambda, "lambda", paste(config.database, temp_output_name))
	dev.off()
}) -> null.dev

dbDisconnect(connection.mysql) -> null.dev;
