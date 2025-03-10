#' Example Dataset: sonar_mines
#'
#' we evaluate the proposed tests through an analysis of a classification task aimed at
#' distinguishing between sonar signals bounced off a metal cylinder and those bounced
#' off a roughly cylindrical rock. The dataset is available at
#' \url{https://archive.ics.uci.edu/dataset/151/connectionist+bench+sonar+mines+vs+rocks}.
#'
#' @format A data frame with 208 rows and 61 variables:
#' \describe{
#'   \item{V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16,V17,V18,V19,V20,V21,V22,V23,
#'   V24,V25,V26,V27,V28,V29,V30,V31,V32,V33,V34,V35,V36,V37,V38,V39,V40,V41,V42,V43,V44,V45,
#'   V46,V47,V48,V49,V50,V51,V52,V53,V54,V55,V56,V57,V58,V59,V60}{Numeric sonar signal attributes (frequencies).}
#'   \item{y}{Class label, a factor with levels 'Mine' and 'Rock'.}
#' }
#' @source from \url{https://archive.ics.uci.edu/dataset/151/connectionist+bench+sonar+mines+vs+rocks}
"sonar_mines"
