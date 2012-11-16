\name{funcXY_J}
\alias{funcXY_J}
\title{Sorts data into array with correct dimensions (JHM)}
\description{
For JHM.
For use with Colonyzer output.
Sorts data into array with correct dimensions.
}
\usage{
funcXY_J(data,data_b,Ma,Mb,N,NoTime_vec,NoSum_vec,NoTime_vec_b,NoSum_vec_b,dimr,dimc)
}
\arguments{
  \item{data}{Growth or Time for control}
  \item{data}{Growth or Time for query}
  \item{M}{M max length of Growth or Time for control}
  \item{M}{M max length of Growth or Time for query}
  \item{N}{Number of ORFs}
  \item{NoTime_vec}{NoTime for control}
  \item{NoSum_vec}{NoSum for control}
  \item{NoTime_vec_b}{NoTime for query}
  \item{NoSum_vec_b}{NoSum for query}
  \item{dimr}{Number of rows in array}
  \item{dimc}{Number of columns in array}
}
\value{
Array
}
\keyword{qfaBayes}

