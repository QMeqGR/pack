def setup(app):
  app.add_crossref_type(
    directivename = "confval",
    rolename      = "confval",
    indextemplate = "pair: %s; Input value",
  )