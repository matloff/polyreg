\name{quizDocs}
\alias{quizDocs}
\alias{quizzes}

\title{
Course quiz documents
}

\description{

This data is suitable for NLP analysis.  It consists of all the quizzes
I've given in undergraduate courses, 143 quizzes in all.  

It is available in two forms.  First, \code{quizzes} is a data.frame,
143 rows and 2 columns.  Row i consists of a single character vector
comprising the entire quiz i, followed by the course name (as an R
factor).  The second form is an R list, 143 elements.  Each list element
is a character vector, one vector element per line of the quiz.

The original documents were LaTeX files.  They have been run through the
\code{detex} utility to remove most LaTeX commands, as well as removing
the LaTeX preambles separately.

The names of the list elements are the course names, as follows:

ECS 50:  a course in machine organization

ECS 132:  an undergraduate course in probabilistic modeling

ECS 145:  a course in scripting languages (Python, R)

ECS 158:  an undergraduate course in parallel computation

ECS 256:  a graduate course in probabilistic modeling

}

