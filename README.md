# Website of Matteo Foglieni

My personal website, running at <https://foglienimatteo.com>. Built with [Quarto](https://quarto.org).

## Preview locally

Install Quarto (<https://quarto.org/docs/get-started/>, or `brew install --cask quarto`), then from this folder:

```bash
quarto preview
```

This opens the site in the browser and reloads it on every save. To build the whole site into `_site/` once:

```bash
quarto render
```

## Write a post

Create `blog/<name-of-the-post>/index.qmd`; the folder name is the address (`/blog/<name-of-the-post>/`). Images go in the same folder.

```yaml
---
title: "Title of the post"
date: 2026-11-13
categories: [math, physics]
---
```

`styleguide/index.qmd` (served at `/styleguide/`, listed nowhere) shows every building block with its source: math, code, executed code, interactive plots, figures, tables, callouts.

## Code that runs

A post may contain `{python}` or `{julia}` cells. They run on your machine when you render, and their output is stored in `_freeze/`, which is committed. The published site is built from those stored results, so the server needs neither Python nor Julia. Only the style guide uses this so far; re-rendering it after an edit needs Python with `jupyter`, `numpy` and `matplotlib`.

`{ojs}` cells (the interactive plots) are JavaScript and run in the reader's browser; they need nothing installed.

## Where things are

| Path | What |
|---|---|
| `_quarto.yml` | site configuration: navigation, footer, page options |
| `index.qmd` | home page |
| `blog/index.qmd` | list of all posts |
| `blog/*/index.qmd` | the posts |
| `Publications/`, `CV/`, `Recommendations/` | the other pages |
| `assets/theme/site.scss` | the look: fonts, colours, spacing |
| `assets/theme/dark.scss` | colours of the dark mode |
| `assets/head.html` | MathJax configuration and macros |
| `assets/post-list.ejs` | template of the date + title post list |
| `assets/pdf-embed.lua` | the `{{< pdf-embed >}}` shortcode (CV, posters) |
| `.github/workflows/publish.yml` | builds and deploys the site on every push to the main branch |

## Deploy

GitHub Actions renders the site and publishes it to GitHub Pages on every push to `main` (or `master`). This needs, once, in the repository settings: **Pages → Build and deployment → Source: GitHub Actions**.
