% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trajectory_inference.R
\name{reverse_trajectory}
\alias{reverse_trajectory}
\title{Reverse a trajectory}
\usage{
reverse_trajectory(trajectory)
}
\arguments{
\item{trajectory}{A trajectory as returned by \code{\link{infer_trajectory}}.}
}
\value{
The same trajectory, but in the other direction.
}
\description{
Since the direction of the trajectory is not specified, the ordering of a trajectory may be inverted using \code{reverse_trajectory}.
}
\examples{
## Generate an example dataset and infer a trajectory through it
dataset <- generate_dataset(num_genes = 200, num_samples = 400, num_groups = 4)
group_name <- dataset$sample_info$group_name
space <- reduce_dimensionality(dataset$expression, ndim = 2)
traj <- infer_trajectory(space)

## Visualise the trajectory
draw_trajectory_plot(space, group_name, path = traj$path)

## Reverse the trajectory
reverse_traj <- reverse_trajectory(traj)
draw_trajectory_plot(space, group_name, path = reverse_traj$path)

plot(traj$time, reverse_traj$time, type = "l")
}
\seealso{
\code{\link{infer_trajectory}}
}
