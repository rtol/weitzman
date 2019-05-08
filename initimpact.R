if (input$radio == "1") {
  imp1 <- 0
  imp2 <- -0.19
  imp6 <- 0
}
else if (input$radio == "2"){
  imp1 <- -0.12
  imp2 <- -0.16
  imp6 <- 0
}
else if (input$radio == "3"){
  imp1 <- 0
  imp2 <- -0.22
  imp6 <- -0.0000371
}
else {
  imp1 <- 0
  imp2 <- 0
  imp6 <- 0
}
