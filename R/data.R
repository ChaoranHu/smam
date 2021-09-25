##' GPS data of f109
##'
##' A dataset of GPS coordinates of a mature female mountain
##' lion living in the Gros Ventre Mountain Range near
##' Jackson, Wyoming. The data were collected by a code-only
##' GPS wildlife tracking collar from 2009 to 2012.
##'
##'
##' @format A data frame with 3919 rows and 4 variables:
##' \describe{
##'     \item{date}{Date when the GPS coordinates were collected}
##'     \item{cumTime}{Standardized time when the GPS coordinates were collected (unit: hour)}
##'     \item{centerE}{Standardized UTM easting (unit: km)}
##'     \item{centerN}{Standardized UTM northing (unit: km)}
##' }
"f109"

##' GPS data of f109 (raw format)
##'
##' The original format of `f109` dataset.
##'
##' @format A data frame with 3917 rows and 3 variables:
##' \describe{
##'     \item{t}{Date and time when the GPS coordinates were collected (unit: year)}
##'     \item{dE}{Standardized UTM easting (unit: meter)}
##'     \item{dN}{Standardized UTM northing (unit: meter)}
##' }
"f109raw"