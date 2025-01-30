library(bulkAnalyseR)
library(ggplot2)

expr.matrix = 'expr_matrix/expr-matrix.csv'
metadata = 'metadata.csv'

expr.matrix = read.csv(file = expr.matrix, row.names = 1)
write.table(expr.matrix, 'processed/raw_expression.csv')

metadata = read.csv(file = metadata)

output.dir = 'bulkanalyser_app_filtered'
app.title = "BulkAnalyseR App"
organism = "mmusculus"
organism.db = "org.Mm.eg.db"

for(i in 2:ncol(metadata)) {
    metadata[,i] <- as.factor(metadata[, i])  
}

expr.matrix = preprocessExpressionMatrix(expr.matrix, output.plot = TRUE)
print(dim(expr.matrix))
expr.matrix <- expr.matrix[apply(expr.matrix, 1, max) > 100, ]
print(dim(expr.matrix))
ggsave('processed/plot.png')
write.table(expr.matrix, 'processed/normalised_expression.csv')


generateShinyApp(
  shiny.dir = output.dir,
  app.title = app.title,
  modality = "RNA",
  expression.matrix = expr.matrix,
  metadata = metadata,
  organism = organism,
  org.db = organism.db 
)
