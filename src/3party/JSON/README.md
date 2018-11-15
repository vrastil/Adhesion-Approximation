





<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
  <link rel="dns-prefetch" href="https://assets-cdn.github.com">
  <link rel="dns-prefetch" href="https://avatars0.githubusercontent.com">
  <link rel="dns-prefetch" href="https://avatars1.githubusercontent.com">
  <link rel="dns-prefetch" href="https://avatars2.githubusercontent.com">
  <link rel="dns-prefetch" href="https://avatars3.githubusercontent.com">
  <link rel="dns-prefetch" href="https://github-cloud.s3.amazonaws.com">
  <link rel="dns-prefetch" href="https://user-images.githubusercontent.com/">



  <link crossorigin="anonymous" media="all" integrity="sha512-vEiCH42J0Z75K8k43QJVo3TeWXpPSjhuyaJLAGYvTHNBexT4HQb1wFP7XY2EKbK37eFbQEk2Z2znKPIUbdMJoA==" rel="stylesheet" href="https://assets-cdn.github.com/assets/frameworks-a92bcfa8d646a4f8874998ac7b7ec3b8.css" />
  <link crossorigin="anonymous" media="all" integrity="sha512-0QmNrnxV5kanLB5PDJ5qKBQZxT6mGBU4ybrofSilTcs5oiQArMz1whpZAyZPdS6aEHK2Sldk1ETxVsajS7Ivbg==" rel="stylesheet" href="https://assets-cdn.github.com/assets/github-9477e1913356a43e26fea78ac0a55d32.css" />
  
  
  <link crossorigin="anonymous" media="all" integrity="sha512-QbxBZqAu/zl5DDhd2QizVMJEVzg31c5Mlqi5uPXG9U7GQFnAk64qWU/ngS64iYhFycI93jHr5hLoAIWvTswe4w==" rel="stylesheet" href="https://assets-cdn.github.com/assets/site-d5d073510fa15e2c00f1b636632a01ff.css" />
  
  

  <meta name="viewport" content="width=device-width">
  
  <title>json/README.md at develop · nlohmann/json · GitHub</title>
    <meta name="description" content="JSON for Modern C++. Contribute to nlohmann/json development by creating an account on GitHub.">
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub">
  <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub">
  <meta property="fb:app_id" content="1401488693436528">

    
    <meta property="og:image" content="https://avatars3.githubusercontent.com/u/159488?s=400&amp;v=4" /><meta property="og:site_name" content="GitHub" /><meta property="og:type" content="object" /><meta property="og:title" content="nlohmann/json" /><meta property="og:url" content="https://github.com/nlohmann/json" /><meta property="og:description" content="JSON for Modern C++. Contribute to nlohmann/json development by creating an account on GitHub." />

  <link rel="assets" href="https://assets-cdn.github.com/">
  
  <meta name="pjax-timeout" content="1000">
  
  <meta name="request-id" content="BB0E:78F0:D39F42F:13032429:5BED5335" data-pjax-transient>


  

  <meta name="selected-link" value="repo_source" data-pjax-transient>

      <meta name="google-site-verification" content="KT5gs8h0wvaagLKAVWq8bbeNwnZZK1r1XQysX3xurLU">
    <meta name="google-site-verification" content="ZzhVyEFwb7w3e0-uOTltm8Jsck2F5StVihD0exw2fsA">
    <meta name="google-site-verification" content="GXs5KoUUkNCoaAZn7wPN-t01Pywp9M3sEjnt_3_ZWPc">

  <meta name="octolytics-host" content="collector.githubapp.com" /><meta name="octolytics-app-id" content="github" /><meta name="octolytics-event-url" content="https://collector.githubapp.com/github-external/browser_event" /><meta name="octolytics-dimension-request_id" content="BB0E:78F0:D39F42F:13032429:5BED5335" /><meta name="octolytics-dimension-region_edge" content="ams" /><meta name="octolytics-dimension-region_render" content="iad" />
<meta name="analytics-location" content="/&lt;user-name&gt;/&lt;repo-name&gt;/blob/show" data-pjax-transient="true" />



    <meta name="google-analytics" content="UA-3769691-2">


<meta class="js-ga-set" name="dimension1" content="Logged Out">



  

      <meta name="hostname" content="github.com">
    <meta name="user-login" content="">

      <meta name="expected-hostname" content="github.com">
    <meta name="js-proxy-site-detection-payload" content="M2E3MmI0MjQ1NTUzODcxMGJlODIxYjBiMzYyYjcxNWNjOGM1ZmFjZGVkMTIwZDlkMDcxNzM5YzA4ZTBjYTZiZnx7InJlbW90ZV9hZGRyZXNzIjoiMTQ3LjIzMS4xOS4zNCIsInJlcXVlc3RfaWQiOiJCQjBFOjc4RjA6RDM5RjQyRjoxMzAzMjQyOTo1QkVENTMzNSIsInRpbWVzdGFtcCI6MTU0MjI3OTk5MCwiaG9zdCI6ImdpdGh1Yi5jb20ifQ==">

    <meta name="enabled-features" content="DASHBOARD_V2_LAYOUT_OPT_IN,EXPLORE_DISCOVER_REPOSITORIES,UNIVERSE_BANNER,MARKETPLACE_PLAN_RESTRICTION_EDITOR">

  <meta name="html-safe-nonce" content="41cfb041bf989dc715ed9862e85cefc55294a097">

  <meta http-equiv="x-pjax-version" content="2cfbed1aa5b3591dd14062561b2c8961">
  

      <link href="https://github.com/nlohmann/json/commits/develop.atom" rel="alternate" title="Recent Commits to json:develop" type="application/atom+xml">

  <meta name="go-import" content="github.com/nlohmann/json git https://github.com/nlohmann/json.git">

  <meta name="octolytics-dimension-user_id" content="159488" /><meta name="octolytics-dimension-user_login" content="nlohmann" /><meta name="octolytics-dimension-repository_id" content="11171548" /><meta name="octolytics-dimension-repository_nwo" content="nlohmann/json" /><meta name="octolytics-dimension-repository_public" content="true" /><meta name="octolytics-dimension-repository_is_fork" content="false" /><meta name="octolytics-dimension-repository_network_root_id" content="11171548" /><meta name="octolytics-dimension-repository_network_root_nwo" content="nlohmann/json" /><meta name="octolytics-dimension-repository_explore_github_marketplace_ci_cta_shown" content="false" />


    <link rel="canonical" href="https://github.com/nlohmann/json/blob/develop/README.md" data-pjax-transient>


  <meta name="browser-stats-url" content="https://api.github.com/_private/browser/stats">

  <meta name="browser-errors-url" content="https://api.github.com/_private/browser/errors">

  <link rel="mask-icon" href="https://assets-cdn.github.com/pinned-octocat.svg" color="#000000">
  <link rel="icon" type="image/x-icon" class="js-site-favicon" href="https://assets-cdn.github.com/favicon.ico">

<meta name="theme-color" content="#1e2327">



  <link rel="manifest" href="/manifest.json" crossOrigin="use-credentials">

  </head>

  <body class="logged-out env-production page-blob">
    

  <div class="position-relative js-header-wrapper ">
    <a href="#start-of-content" tabindex="1" class="px-2 py-4 bg-blue text-white show-on-focus js-skip-to-content">Skip to content</a>
    <div id="js-pjax-loader-bar" class="pjax-loader-bar"><div class="progress"></div></div>

    
    
    



        
<header class="Header header-logged-out  position-relative f4 py-3" role="banner">
  <div class="container-lg d-flex px-3">
    <div class="d-flex flex-justify-between flex-items-center">
      <a class="header-logo-invertocat my-0" href="https://github.com/" aria-label="Homepage" data-ga-click="(Logged out) Header, go to homepage, icon:logo-wordmark">
        <svg height="32" class="octicon octicon-mark-github" viewBox="0 0 16 16" version="1.1" width="32" aria-hidden="true"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"/></svg>
      </a>

    </div>

    <div class="HeaderMenu d-flex flex-justify-between flex-auto">
        <nav class="mt-0">
          <ul class="d-flex list-style-none">
              <li class="ml-2">
                <a class="js-selected-navigation-item HeaderNavlink px-0 py-2 m-0" data-ga-click="Header, click, Nav menu - item:features" data-selected-links="/features /features/project-management /features/code-review /features/project-management /features/integrations /features" href="/features">
                  Features
</a>              </li>
              <li class="ml-4">
                <a class="js-selected-navigation-item HeaderNavlink px-0 py-2 m-0" data-ga-click="Header, click, Nav menu - item:business" data-selected-links="/business /business/security /business/customers /business" href="/business">
                  Business
</a>              </li>

              <li class="ml-4">
                <a class="js-selected-navigation-item HeaderNavlink px-0 py-2 m-0" data-ga-click="Header, click, Nav menu - item:explore" data-selected-links="/explore /trending /trending/developers /integrations /integrations/feature/code /integrations/feature/collaborate /integrations/feature/ship showcases showcases_search showcases_landing /explore" href="/explore">
                  Explore
</a>              </li>

              <li class="ml-4">
                    <a class="js-selected-navigation-item HeaderNavlink px-0 py-2 m-0" data-ga-click="Header, click, Nav menu - item:marketplace" data-selected-links=" /marketplace" href="/marketplace">
                      Marketplace
</a>              </li>
              <li class="ml-4">
                <a class="js-selected-navigation-item HeaderNavlink px-0 py-2 m-0" data-ga-click="Header, click, Nav menu - item:pricing" data-selected-links="/pricing /pricing/developer /pricing/team /pricing/business-hosted /pricing/business-enterprise /pricing" href="/pricing">
                  Pricing
</a>              </li>
          </ul>
        </nav>

      <div class="d-flex">
          <div class="d-lg-flex flex-items-center mr-3">
            <div class="header-search scoped-search site-scoped-search js-site-search position-relative js-jump-to"
  role="combobox"
  aria-owns="jump-to-results"
  aria-label="Search or jump to"
  aria-haspopup="listbox"
  aria-expanded="false"
>
  <div class="position-relative">
    <!-- '"` --><!-- </textarea></xmp> --></option></form><form class="js-site-search-form" data-scope-type="Repository" data-scope-id="11171548" data-scoped-search-url="/nlohmann/json/search" data-unscoped-search-url="/search" action="/nlohmann/json/search" accept-charset="UTF-8" method="get"><input name="utf8" type="hidden" value="&#x2713;" />
      <label class="form-control header-search-wrapper header-search-wrapper-jump-to position-relative d-flex flex-justify-between flex-items-center js-chromeless-input-container">
        <input type="text"
          class="form-control header-search-input jump-to-field js-jump-to-field js-site-search-focus js-site-search-field is-clearable"
          data-hotkey="s,/"
          name="q"
          value=""
          placeholder="Search"
          data-unscoped-placeholder="Search GitHub"
          data-scoped-placeholder="Search"
          autocapitalize="off"
          aria-autocomplete="list"
          aria-controls="jump-to-results"
          aria-label="Search"
          data-jump-to-suggestions-path="/_graphql/GetSuggestedNavigationDestinations#csrf-token=W29coikEF10WWimN+2zURT/zKjau79vR2H2vI+YbbgNFVLkmCC2tThhqSJ/6DwEvq/OuoyXgAnbRC5vjUjHXow=="
          spellcheck="false"
          autocomplete="off"
          >
          <input type="hidden" class="js-site-search-type-field" name="type" >
            <img src="https://assets-cdn.github.com/images/search-shortcut-hint.svg" alt="" class="mr-2 header-search-key-slash">

            <div class="Box position-absolute overflow-hidden d-none jump-to-suggestions js-jump-to-suggestions-container">
              <ul class="d-none js-jump-to-suggestions-template-container">
                <li class="d-flex flex-justify-start flex-items-center p-0 f5 navigation-item js-navigation-item" role="option">
                  <a tabindex="-1" class="no-underline d-flex flex-auto flex-items-center p-2 jump-to-suggestions-path js-jump-to-suggestion-path js-navigation-open" href="">
                    <div class="jump-to-octicon js-jump-to-octicon flex-shrink-0 mr-2 text-center d-none">
                      <svg height="16" width="16" class="octicon octicon-repo flex-shrink-0 js-jump-to-octicon-repo d-none" title="Repository" aria-label="Repository" viewBox="0 0 12 16" version="1.1" role="img"><path fill-rule="evenodd" d="M4 9H3V8h1v1zm0-3H3v1h1V6zm0-2H3v1h1V4zm0-2H3v1h1V2zm8-1v12c0 .55-.45 1-1 1H6v2l-1.5-1.5L3 16v-2H1c-.55 0-1-.45-1-1V1c0-.55.45-1 1-1h10c.55 0 1 .45 1 1zm-1 10H1v2h2v-1h3v1h5v-2zm0-10H2v9h9V1z"/></svg>
                      <svg height="16" width="16" class="octicon octicon-project flex-shrink-0 js-jump-to-octicon-project d-none" title="Project" aria-label="Project" viewBox="0 0 15 16" version="1.1" role="img"><path fill-rule="evenodd" d="M10 12h3V2h-3v10zm-4-2h3V2H6v8zm-4 4h3V2H2v12zm-1 1h13V1H1v14zM14 0H1a1 1 0 0 0-1 1v14a1 1 0 0 0 1 1h13a1 1 0 0 0 1-1V1a1 1 0 0 0-1-1z"/></svg>
                      <svg height="16" width="16" class="octicon octicon-search flex-shrink-0 js-jump-to-octicon-search d-none" title="Search" aria-label="Search" viewBox="0 0 16 16" version="1.1" role="img"><path fill-rule="evenodd" d="M15.7 13.3l-3.81-3.83A5.93 5.93 0 0 0 13 6c0-3.31-2.69-6-6-6S1 2.69 1 6s2.69 6 6 6c1.3 0 2.48-.41 3.47-1.11l3.83 3.81c.19.2.45.3.7.3.25 0 .52-.09.7-.3a.996.996 0 0 0 0-1.41v.01zM7 10.7c-2.59 0-4.7-2.11-4.7-4.7 0-2.59 2.11-4.7 4.7-4.7 2.59 0 4.7 2.11 4.7 4.7 0 2.59-2.11 4.7-4.7 4.7z"/></svg>
                    </div>

                    <img class="avatar mr-2 flex-shrink-0 js-jump-to-suggestion-avatar d-none" alt="" aria-label="Team" src="" width="28" height="28">

                    <div class="jump-to-suggestion-name js-jump-to-suggestion-name flex-auto overflow-hidden text-left no-wrap css-truncate css-truncate-target">
                    </div>

                    <div class="border rounded-1 flex-shrink-0 bg-gray px-1 text-gray-light ml-1 f6 d-none js-jump-to-badge-search">
                      <span class="js-jump-to-badge-search-text-default d-none" aria-label="in this repository">
                        In this repository
                      </span>
                      <span class="js-jump-to-badge-search-text-global d-none" aria-label="in all of GitHub">
                        All GitHub
                      </span>
                      <span aria-hidden="true" class="d-inline-block ml-1 v-align-middle">↵</span>
                    </div>

                    <div aria-hidden="true" class="border rounded-1 flex-shrink-0 bg-gray px-1 text-gray-light ml-1 f6 d-none d-on-nav-focus js-jump-to-badge-jump">
                      Jump to
                      <span class="d-inline-block ml-1 v-align-middle">↵</span>
                    </div>
                  </a>
                </li>
              </ul>
              <ul class="d-none js-jump-to-no-results-template-container">
                <li class="d-flex flex-justify-center flex-items-center p-3 f5 d-none">
                  <span class="text-gray">No suggested jump to results</span>
                </li>
              </ul>

              <ul id="jump-to-results" role="listbox" class="js-navigation-container jump-to-suggestions-results-container js-jump-to-suggestions-results-container" >
                <li class="d-flex flex-justify-center flex-items-center p-0 f5">
                  <img src="https://assets-cdn.github.com/images/spinners/octocat-spinner-128.gif" alt="Octocat Spinner Icon" class="m-2" width="28">
                </li>
              </ul>
            </div>
      </label>
</form>  </div>
</div>

          </div>

        <span class="d-inline-block">
            <div class="HeaderNavlink px-0 py-2 m-0">
              <a class="text-bold text-white no-underline" href="/login?return_to=%2Fnlohmann%2Fjson%2Fblob%2Fdevelop%2FREADME.md" data-ga-click="(Logged out) Header, clicked Sign in, text:sign-in">Sign in</a>
                <span class="text-gray">or</span>
                <a class="text-bold text-white no-underline" href="/join?source=header-repo" data-ga-click="(Logged out) Header, clicked Sign up, text:sign-up">Sign up</a>
            </div>
        </span>
      </div>
    </div>
  </div>
</header>

  </div>

  <div id="start-of-content" class="show-on-focus"></div>

    <div id="js-flash-container">


</div>



  <div role="main" class="application-main ">
        <div itemscope itemtype="http://schema.org/SoftwareSourceCode" class="">
    <div id="js-repo-pjax-container" data-pjax-container >
      







  <div class="pagehead repohead instapaper_ignore readability-menu experiment-repo-nav  ">
    <div class="repohead-details-container clearfix container">

      <ul class="pagehead-actions">
  <li>
      <a href="/login?return_to=%2Fnlohmann%2Fjson"
    class="btn btn-sm btn-with-count tooltipped tooltipped-s"
    aria-label="You must be signed in to watch a repository" rel="nofollow">
    <svg class="octicon octicon-eye v-align-text-bottom" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M8.06 2C3 2 0 8 0 8s3 6 8.06 6C13 14 16 8 16 8s-3-6-7.94-6zM8 12c-2.2 0-4-1.78-4-4 0-2.2 1.8-4 4-4 2.22 0 4 1.8 4 4 0 2.22-1.78 4-4 4zm2-4c0 1.11-.89 2-2 2-1.11 0-2-.89-2-2 0-1.11.89-2 2-2 1.11 0 2 .89 2 2z"/></svg>
    Watch
  </a>
  <a class="social-count" href="/nlohmann/json/watchers"
     aria-label="559 users are watching this repository">
    559
  </a>

  </li>

  <li>
      <a href="/login?return_to=%2Fnlohmann%2Fjson"
    class="btn btn-sm btn-with-count tooltipped tooltipped-s"
    aria-label="You must be signed in to star a repository" rel="nofollow">
    <svg class="octicon octicon-star v-align-text-bottom" viewBox="0 0 14 16" version="1.1" width="14" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M14 6l-4.9-.64L7 1 4.9 5.36 0 6l3.6 3.26L2.67 14 7 11.67 11.33 14l-.93-4.74L14 6z"/></svg>
    Star
  </a>

    <a class="social-count js-social-count" href="/nlohmann/json/stargazers"
      aria-label="11765 users starred this repository">
      11,765
    </a>

  </li>

  <li>
      <a href="/login?return_to=%2Fnlohmann%2Fjson"
        class="btn btn-sm btn-with-count tooltipped tooltipped-s"
        aria-label="You must be signed in to fork a repository" rel="nofollow">
        <svg class="octicon octicon-repo-forked v-align-text-bottom" viewBox="0 0 10 16" version="1.1" width="10" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M8 1a1.993 1.993 0 0 0-1 3.72V6L5 8 3 6V4.72A1.993 1.993 0 0 0 2 1a1.993 1.993 0 0 0-1 3.72V6.5l3 3v1.78A1.993 1.993 0 0 0 5 15a1.993 1.993 0 0 0 1-3.72V9.5l3-3V4.72A1.993 1.993 0 0 0 8 1zM2 4.2C1.34 4.2.8 3.65.8 3c0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zm3 10c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zm3-10c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2z"/></svg>
        Fork
      </a>

    <a href="/nlohmann/json/network/members" class="social-count"
       aria-label="1892 users forked this repository">
      1,892
    </a>
  </li>
</ul>

      <h1 class="public ">
  <svg class="octicon octicon-repo" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9H3V8h1v1zm0-3H3v1h1V6zm0-2H3v1h1V4zm0-2H3v1h1V2zm8-1v12c0 .55-.45 1-1 1H6v2l-1.5-1.5L3 16v-2H1c-.55 0-1-.45-1-1V1c0-.55.45-1 1-1h10c.55 0 1 .45 1 1zm-1 10H1v2h2v-1h3v1h5v-2zm0-10H2v9h9V1z"/></svg>
  <span class="author" itemprop="author"><a class="url fn" rel="author" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=159488" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann">nlohmann</a></span><!--
--><span class="path-divider">/</span><!--
--><strong itemprop="name"><a data-pjax="#js-repo-pjax-container" href="/nlohmann/json">json</a></strong>

</h1>

    </div>
    
<nav class="reponav js-repo-nav js-sidenav-container-pjax container"
     itemscope
     itemtype="http://schema.org/BreadcrumbList"
     role="navigation"
     data-pjax="#js-repo-pjax-container">

  <span itemscope itemtype="http://schema.org/ListItem" itemprop="itemListElement">
    <a class="js-selected-navigation-item selected reponav-item" itemprop="url" data-hotkey="g c" aria-current="page" data-selected-links="repo_source repo_downloads repo_commits repo_releases repo_tags repo_branches repo_packages /nlohmann/json" href="/nlohmann/json">
      <svg class="octicon octicon-code" viewBox="0 0 14 16" version="1.1" width="14" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M9.5 3L8 4.5 11.5 8 8 11.5 9.5 13 14 8 9.5 3zm-5 0L0 8l4.5 5L6 11.5 2.5 8 6 4.5 4.5 3z"/></svg>
      <span itemprop="name">Code</span>
      <meta itemprop="position" content="1">
</a>  </span>

    <span itemscope itemtype="http://schema.org/ListItem" itemprop="itemListElement">
      <a itemprop="url" data-hotkey="g i" class="js-selected-navigation-item reponav-item" data-selected-links="repo_issues repo_labels repo_milestones /nlohmann/json/issues" href="/nlohmann/json/issues">
        <svg class="octicon octicon-issue-opened" viewBox="0 0 14 16" version="1.1" width="14" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7 2.3c3.14 0 5.7 2.56 5.7 5.7s-2.56 5.7-5.7 5.7A5.71 5.71 0 0 1 1.3 8c0-3.14 2.56-5.7 5.7-5.7zM7 1C3.14 1 0 4.14 0 8s3.14 7 7 7 7-3.14 7-7-3.14-7-7-7zm1 3H6v5h2V4zm0 6H6v2h2v-2z"/></svg>
        <span itemprop="name">Issues</span>
        <span class="Counter">21</span>
        <meta itemprop="position" content="2">
</a>    </span>

  <span itemscope itemtype="http://schema.org/ListItem" itemprop="itemListElement">
    <a data-hotkey="g p" itemprop="url" class="js-selected-navigation-item reponav-item" data-selected-links="repo_pulls checks /nlohmann/json/pulls" href="/nlohmann/json/pulls">
      <svg class="octicon octicon-git-pull-request" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M11 11.28V5c-.03-.78-.34-1.47-.94-2.06C9.46 2.35 8.78 2.03 8 2H7V0L4 3l3 3V4h1c.27.02.48.11.69.31.21.2.3.42.31.69v6.28A1.993 1.993 0 0 0 10 15a1.993 1.993 0 0 0 1-3.72zm-1 2.92c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zM4 3c0-1.11-.89-2-2-2a1.993 1.993 0 0 0-1 3.72v6.56A1.993 1.993 0 0 0 2 15a1.993 1.993 0 0 0 1-3.72V4.72c.59-.34 1-.98 1-1.72zm-.8 10c0 .66-.55 1.2-1.2 1.2-.65 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2zM2 4.2C1.34 4.2.8 3.65.8 3c0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2z"/></svg>
      <span itemprop="name">Pull requests</span>
      <span class="Counter">5</span>
      <meta itemprop="position" content="3">
</a>  </span>





  <a class="js-selected-navigation-item reponav-item" data-selected-links="repo_graphs repo_contributors dependency_graph pulse alerts /nlohmann/json/pulse" href="/nlohmann/json/pulse">
    <svg class="octicon octicon-graph" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M16 14v1H0V0h1v14h15zM5 13H3V8h2v5zm4 0H7V3h2v10zm4 0h-2V6h2v7z"/></svg>
    Insights
</a>

</nav>


  </div>

<div class="container new-discussion-timeline experiment-repo-nav  ">
  <div class="repository-content ">

    

  
    <a class="d-none js-permalink-shortcut" data-hotkey="y" href="/nlohmann/json/blob/da81e7be22d431dafaf7618611a5513d71f5fafc/README.md">Permalink</a>

    <!-- blob contrib key: blob_contributors:v21:90f82429840a1f2252970b290543195a -->

        <div class="signup-prompt-bg rounded-1">
      <div class="signup-prompt p-4 text-center mb-4 rounded-1">
        <div class="position-relative">
          <!-- '"` --><!-- </textarea></xmp> --></option></form><form action="/site/dismiss_signup_prompt" accept-charset="UTF-8" method="post"><input name="utf8" type="hidden" value="&#x2713;" /><input type="hidden" name="authenticity_token" value="DDhfyfWZQKLnvw2ncQOs1M3yu/Uv0eFqLoUlv8LSfMzJoxDJwvx+29Jt1IteBBFA8taQN4o7qiHeipuKyS/msw==" />
            <button type="submit" class="position-absolute top-0 right-0 btn-link link-gray" data-ga-click="(Logged out) Sign up prompt, clicked Dismiss, text:dismiss">
              Dismiss
            </button>
</form>          <h3 class="pt-2">Join GitHub today</h3>
          <p class="col-6 mx-auto">GitHub is home to over 28 million developers working together to host and review code, manage projects, and build software together.</p>
          <a class="btn btn-primary" href="/join?source=prompt-blob-show" data-ga-click="(Logged out) Sign up prompt, clicked Sign up, text:sign-up">Sign up</a>
        </div>
      </div>
    </div>


    <div class="file-navigation">
      
<div class="select-menu branch-select-menu js-menu-container js-select-menu float-left">
  <button class=" btn btn-sm select-menu-button js-menu-target css-truncate" data-hotkey="w"
    
    type="button" aria-label="Switch branches or tags" aria-expanded="false" aria-haspopup="true">
      <i>Branch:</i>
      <span class="js-select-button css-truncate-target">develop</span>
  </button>

  <div class="select-menu-modal-holder js-menu-content js-navigation-container" data-pjax>

    <div class="select-menu-modal">
      <div class="select-menu-header">
        <svg class="octicon octicon-x js-menu-close" role="img" aria-label="Close" viewBox="0 0 12 16" version="1.1" width="12" height="16"><path fill-rule="evenodd" d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48L7.48 8z"/></svg>
        <span class="select-menu-title">Switch branches/tags</span>
      </div>

      <div class="select-menu-filters">
        <div class="select-menu-text-filter">
          <input type="text" aria-label="Filter branches/tags" id="context-commitish-filter-field" class="form-control js-filterable-field js-navigation-enable" placeholder="Filter branches/tags">
        </div>
        <div class="select-menu-tabs" role="tablist">
          <ul>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="branches" data-filter-placeholder="Filter branches/tags" class="js-select-menu-tab" role="tab">Branches</a>
            </li>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="tags" data-filter-placeholder="Find a tag…" class="js-select-menu-tab" role="tab">Tags</a>
            </li>
          </ul>
        </div>
      </div>

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="branches" role="menu">

        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/nlohmann/json/blob/coverity_scan/README.md"
               data-name="coverity_scan"
               data-skip-pjax="true"
               rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                coverity_scan
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open selected"
               href="/nlohmann/json/blob/develop/README.md"
               data-name="develop"
               data-skip-pjax="true"
               rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                develop
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/nlohmann/json/blob/gh-pages/README.md"
               data-name="gh-pages"
               data-skip-pjax="true"
               rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                gh-pages
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
               href="/nlohmann/json/blob/master/README.md"
               data-name="master"
               data-skip-pjax="true"
               rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                master
              </span>
            </a>
        </div>

          <div class="select-menu-no-results">Nothing to show</div>
      </div>

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="tags">
        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v3.4.0/README.md"
              data-name="v3.4.0"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v3.4.0">
                v3.4.0
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v3.3.0/README.md"
              data-name="v3.3.0"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v3.3.0">
                v3.3.0
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v3.2.0/README.md"
              data-name="v3.2.0"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v3.2.0">
                v3.2.0
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v3.1.2/README.md"
              data-name="v3.1.2"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v3.1.2">
                v3.1.2
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v3.1.1/README.md"
              data-name="v3.1.1"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v3.1.1">
                v3.1.1
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v3.1.0/README.md"
              data-name="v3.1.0"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v3.1.0">
                v3.1.0
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v3.0.1/README.md"
              data-name="v3.0.1"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v3.0.1">
                v3.0.1
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v3.0.0/README.md"
              data-name="v3.0.0"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v3.0.0">
                v3.0.0
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v2.1.1/README.md"
              data-name="v2.1.1"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v2.1.1">
                v2.1.1
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v2.1.0/README.md"
              data-name="v2.1.0"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v2.1.0">
                v2.1.0
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v2.0.10/README.md"
              data-name="v2.0.10"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v2.0.10">
                v2.0.10
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v2.0.9/README.md"
              data-name="v2.0.9"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v2.0.9">
                v2.0.9
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v2.0.8/README.md"
              data-name="v2.0.8"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v2.0.8">
                v2.0.8
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v2.0.7/README.md"
              data-name="v2.0.7"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v2.0.7">
                v2.0.7
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v2.0.6/README.md"
              data-name="v2.0.6"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v2.0.6">
                v2.0.6
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v2.0.5/README.md"
              data-name="v2.0.5"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v2.0.5">
                v2.0.5
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v2.0.4/README.md"
              data-name="v2.0.4"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v2.0.4">
                v2.0.4
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v2.0.3/README.md"
              data-name="v2.0.3"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v2.0.3">
                v2.0.3
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v2.0.2/README.md"
              data-name="v2.0.2"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v2.0.2">
                v2.0.2
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v2.0.1/README.md"
              data-name="v2.0.1"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v2.0.1">
                v2.0.1
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v2.0.0/README.md"
              data-name="v2.0.0"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v2.0.0">
                v2.0.0
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v1.1.0/README.md"
              data-name="v1.1.0"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v1.1.0">
                v1.1.0
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v1.0.0/README.md"
              data-name="v1.0.0"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v1.0.0">
                v1.0.0
              </span>
            </a>
            <a class="select-menu-item js-navigation-item js-navigation-open "
              href="/nlohmann/json/tree/v1.0.0-rc1/README.md"
              data-name="v1.0.0-rc1"
              data-skip-pjax="true"
              rel="nofollow">
              <svg class="octicon octicon-check select-menu-item-icon" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5L12 5z"/></svg>
              <span class="select-menu-item-text css-truncate-target" title="v1.0.0-rc1">
                v1.0.0-rc1
              </span>
            </a>
        </div>

        <div class="select-menu-no-results">Nothing to show</div>
      </div>

    </div>
  </div>
</div>

      <div class="BtnGroup float-right">
        <a href="/nlohmann/json/find/develop"
              class="js-pjax-capture-input btn btn-sm BtnGroup-item"
              data-pjax
              data-hotkey="t">
          Find file
        </a>
        <clipboard-copy for="blob-path" class="btn btn-sm BtnGroup-item">
          Copy path
        </clipboard-copy>
      </div>
      <div id="blob-path" class="breadcrumb">
        <span class="repo-root js-repo-root"><span class="js-path-segment"><a data-pjax="true" href="/nlohmann/json"><span>json</span></a></span></span><span class="separator">/</span><strong class="final-path">README.md</strong>
      </div>
    </div>


    
  <div class="commit-tease">
      <span class="float-right">
        <a class="commit-tease-sha" href="/nlohmann/json/commit/689382a722e95decd6a1ec855ce263e42b503f7c" data-pjax>
          689382a
        </a>
        <relative-time datetime="2018-11-02T08:35:17Z">Nov 2, 2018</relative-time>
      </span>
      <div>
        <a rel="contributor" data-skip-pjax="true" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=1353258" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/ax3l"><img class="avatar" src="https://avatars2.githubusercontent.com/u/1353258?s=40&amp;v=4" width="20" height="20" alt="@ax3l" /></a>
        <a class="user-mention" rel="contributor" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=1353258" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/ax3l">ax3l</a>
          <a data-pjax="true" title="Fix EOL Whitespaces &amp; CMake Spelling

Fix little leftover EOL whitespaces in `CMakeLists.txt` and
a spelling of CMake in README.md" class="message" href="/nlohmann/json/commit/689382a722e95decd6a1ec855ce263e42b503f7c">Fix EOL Whitespaces &amp; CMake Spelling</a>
      </div>

    <div class="commit-tease-contributors">
      
<details class="details-reset details-overlay details-overlay-dark lh-default text-gray-dark float-left mr-2" id="blob_contributors_box">
  <summary class="btn-link" aria-haspopup="dialog"  >
    
    <span><strong>38</strong> contributors</span>
  </summary>
  <details-dialog class="Box Box--overlay d-flex flex-column anim-fade-in fast " aria-label="Users who have contributed to this file">
    <div class="Box-header">
      <button class="Box-btn-octicon btn-octicon float-right" type="button" aria-label="Close dialog" data-close-dialog>
        <svg class="octicon octicon-x" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48L7.48 8z"/></svg>
      </button>
      <h3 class="Box-title">Users who have contributed to this file</h3>
    </div>
    
        <ul class="list-style-none overflow-auto">
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/nlohmann">
                <img class="avatar mr-2" alt="" src="https://avatars0.githubusercontent.com/u/159488?s=40&amp;v=4" width="20" height="20" />
                nlohmann
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/pfultz2">
                <img class="avatar mr-2" alt="" src="https://avatars1.githubusercontent.com/u/1306044?s=40&amp;v=4" width="20" height="20" />
                pfultz2
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/patrikhuber">
                <img class="avatar mr-2" alt="" src="https://avatars3.githubusercontent.com/u/4967343?s=40&amp;v=4" width="20" height="20" />
                patrikhuber
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/Type1J">
                <img class="avatar mr-2" alt="" src="https://avatars1.githubusercontent.com/u/413028?s=40&amp;v=4" width="20" height="20" />
                Type1J
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/theodelrieu">
                <img class="avatar mr-2" alt="" src="https://avatars3.githubusercontent.com/u/15652306?s=40&amp;v=4" width="20" height="20" />
                theodelrieu
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/Teemperor">
                <img class="avatar mr-2" alt="" src="https://avatars1.githubusercontent.com/u/736001?s=40&amp;v=4" width="20" height="20" />
                Teemperor
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/jrakow">
                <img class="avatar mr-2" alt="" src="https://avatars2.githubusercontent.com/u/19510812?s=40&amp;v=4" width="20" height="20" />
                jrakow
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/jowr">
                <img class="avatar mr-2" alt="" src="https://avatars1.githubusercontent.com/u/769593?s=40&amp;v=4" width="20" height="20" />
                jowr
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/Itja">
                <img class="avatar mr-2" alt="" src="https://avatars0.githubusercontent.com/u/11618623?s=40&amp;v=4" width="20" height="20" />
                Itja
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/chuckatkins">
                <img class="avatar mr-2" alt="" src="https://avatars1.githubusercontent.com/u/320135?s=40&amp;v=4" width="20" height="20" />
                chuckatkins
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/ax3l">
                <img class="avatar mr-2" alt="" src="https://avatars2.githubusercontent.com/u/1353258?s=40&amp;v=4" width="20" height="20" />
                ax3l
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/wla80">
                <img class="avatar mr-2" alt="" src="https://avatars0.githubusercontent.com/u/13490408?s=40&amp;v=4" width="20" height="20" />
                wla80
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/vasild">
                <img class="avatar mr-2" alt="" src="https://avatars2.githubusercontent.com/u/266751?s=40&amp;v=4" width="20" height="20" />
                vasild
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/vog">
                <img class="avatar mr-2" alt="" src="https://avatars2.githubusercontent.com/u/412749?s=40&amp;v=4" width="20" height="20" />
                vog
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/Pipeliner">
                <img class="avatar mr-2" alt="" src="https://avatars3.githubusercontent.com/u/598225?s=40&amp;v=4" width="20" height="20" />
                Pipeliner
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/feroldi">
                <img class="avatar mr-2" alt="" src="https://avatars2.githubusercontent.com/u/8634526?s=40&amp;v=4" width="20" height="20" />
                feroldi
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/abolz">
                <img class="avatar mr-2" alt="" src="https://avatars1.githubusercontent.com/u/1213085?s=40&amp;v=4" width="20" height="20" />
                abolz
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/koponomarenko">
                <img class="avatar mr-2" alt="" src="https://avatars2.githubusercontent.com/u/3629232?s=40&amp;v=4" width="20" height="20" />
                koponomarenko
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/kevin--">
                <img class="avatar mr-2" alt="" src="https://avatars0.githubusercontent.com/u/3334774?s=40&amp;v=4" width="20" height="20" />
                kevin--
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/jammehcow">
                <img class="avatar mr-2" alt="" src="https://avatars3.githubusercontent.com/u/4962764?s=40&amp;v=4" width="20" height="20" />
                jammehcow
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/jaredgrubb">
                <img class="avatar mr-2" alt="" src="https://avatars3.githubusercontent.com/u/1256381?s=40&amp;v=4" width="20" height="20" />
                jaredgrubb
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/henryiii">
                <img class="avatar mr-2" alt="" src="https://avatars0.githubusercontent.com/u/4616906?s=40&amp;v=4" width="20" height="20" />
                henryiii
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/gregmarr">
                <img class="avatar mr-2" alt="" src="https://avatars0.githubusercontent.com/u/8569738?s=40&amp;v=4" width="20" height="20" />
                gregmarr
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/zerodefect">
                <img class="avatar mr-2" alt="" src="https://avatars2.githubusercontent.com/u/26778249?s=40&amp;v=4" width="20" height="20" />
                zerodefect
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/Dobiasd">
                <img class="avatar mr-2" alt="" src="https://avatars0.githubusercontent.com/u/5544610?s=40&amp;v=4" width="20" height="20" />
                Dobiasd
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/dan-42">
                <img class="avatar mr-2" alt="" src="https://avatars3.githubusercontent.com/u/1706857?s=40&amp;v=4" width="20" height="20" />
                dan-42
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/daixtrose">
                <img class="avatar mr-2" alt="" src="https://avatars2.githubusercontent.com/u/5588692?s=40&amp;v=4" width="20" height="20" />
                daixtrose
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/seeekr">
                <img class="avatar mr-2" alt="" src="https://avatars1.githubusercontent.com/u/302886?s=40&amp;v=4" width="20" height="20" />
                seeekr
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/ChrisKitching">
                <img class="avatar mr-2" alt="" src="https://avatars3.githubusercontent.com/u/3526849?s=40&amp;v=4" width="20" height="20" />
                ChrisKitching
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/berkus">
                <img class="avatar mr-2" alt="" src="https://avatars0.githubusercontent.com/u/80306?s=40&amp;v=4" width="20" height="20" />
                berkus
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/aqnouch">
                <img class="avatar mr-2" alt="" src="https://avatars0.githubusercontent.com/u/3692812?s=40&amp;v=4" width="20" height="20" />
                aqnouch
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/Annihil">
                <img class="avatar mr-2" alt="" src="https://avatars0.githubusercontent.com/u/16704309?s=40&amp;v=4" width="20" height="20" />
                Annihil
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/andoma">
                <img class="avatar mr-2" alt="" src="https://avatars2.githubusercontent.com/u/216384?s=40&amp;v=4" width="20" height="20" />
                andoma
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/alex-weej">
                <img class="avatar mr-2" alt="" src="https://avatars3.githubusercontent.com/u/15476?s=40&amp;v=4" width="20" height="20" />
                alex-weej
</a>            </li>
            <li class="Box-row">
              <a class="link-gray-dark no-underline" href="/martin-mfg">
                <img class="avatar mr-2" alt="" src="https://avatars2.githubusercontent.com/u/2026226?s=40&amp;v=4" width="20" height="20" />
                martin-mfg
</a>            </li>
        </ul>

  </details-dialog>
</details>
          <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=159488" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=nlohmann">
      <img class="avatar" src="https://avatars0.githubusercontent.com/u/159488?s=40&amp;v=4" width="20" height="20" alt="@nlohmann" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=1306044" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=pfultz2">
      <img class="avatar" src="https://avatars1.githubusercontent.com/u/1306044?s=40&amp;v=4" width="20" height="20" alt="@pfultz2" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=4967343" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=patrikhuber">
      <img class="avatar" src="https://avatars3.githubusercontent.com/u/4967343?s=40&amp;v=4" width="20" height="20" alt="@patrikhuber" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=413028" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=Type1J">
      <img class="avatar" src="https://avatars1.githubusercontent.com/u/413028?s=40&amp;v=4" width="20" height="20" alt="@Type1J" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=15652306" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=theodelrieu">
      <img class="avatar" src="https://avatars3.githubusercontent.com/u/15652306?s=40&amp;v=4" width="20" height="20" alt="@theodelrieu" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=736001" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=Teemperor">
      <img class="avatar" src="https://avatars1.githubusercontent.com/u/736001?s=40&amp;v=4" width="20" height="20" alt="@Teemperor" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=19510812" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=jrakow">
      <img class="avatar" src="https://avatars2.githubusercontent.com/u/19510812?s=40&amp;v=4" width="20" height="20" alt="@jrakow" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=769593" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=jowr">
      <img class="avatar" src="https://avatars1.githubusercontent.com/u/769593?s=40&amp;v=4" width="20" height="20" alt="@jowr" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=11618623" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=Itja">
      <img class="avatar" src="https://avatars0.githubusercontent.com/u/11618623?s=40&amp;v=4" width="20" height="20" alt="@Itja" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=320135" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=chuckatkins">
      <img class="avatar" src="https://avatars1.githubusercontent.com/u/320135?s=40&amp;v=4" width="20" height="20" alt="@chuckatkins" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=1353258" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=ax3l">
      <img class="avatar" src="https://avatars2.githubusercontent.com/u/1353258?s=40&amp;v=4" width="20" height="20" alt="@ax3l" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=13490408" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=wla80">
      <img class="avatar" src="https://avatars0.githubusercontent.com/u/13490408?s=40&amp;v=4" width="20" height="20" alt="@wla80" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=266751" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=vasild">
      <img class="avatar" src="https://avatars2.githubusercontent.com/u/266751?s=40&amp;v=4" width="20" height="20" alt="@vasild" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=412749" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=vog">
      <img class="avatar" src="https://avatars2.githubusercontent.com/u/412749?s=40&amp;v=4" width="20" height="20" alt="@vog" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=598225" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=Pipeliner">
      <img class="avatar" src="https://avatars3.githubusercontent.com/u/598225?s=40&amp;v=4" width="20" height="20" alt="@Pipeliner" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=8634526" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=feroldi">
      <img class="avatar" src="https://avatars2.githubusercontent.com/u/8634526?s=40&amp;v=4" width="20" height="20" alt="@feroldi" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=1213085" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=abolz">
      <img class="avatar" src="https://avatars1.githubusercontent.com/u/1213085?s=40&amp;v=4" width="20" height="20" alt="@abolz" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=3629232" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=koponomarenko">
      <img class="avatar" src="https://avatars2.githubusercontent.com/u/3629232?s=40&amp;v=4" width="20" height="20" alt="@koponomarenko" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=3334774" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=kevin--">
      <img class="avatar" src="https://avatars0.githubusercontent.com/u/3334774?s=40&amp;v=4" width="20" height="20" alt="@kevin--" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=4962764" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=jammehcow">
      <img class="avatar" src="https://avatars3.githubusercontent.com/u/4962764?s=40&amp;v=4" width="20" height="20" alt="@jammehcow" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=1256381" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=jaredgrubb">
      <img class="avatar" src="https://avatars3.githubusercontent.com/u/1256381?s=40&amp;v=4" width="20" height="20" alt="@jaredgrubb" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=4616906" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=henryiii">
      <img class="avatar" src="https://avatars0.githubusercontent.com/u/4616906?s=40&amp;v=4" width="20" height="20" alt="@henryiii" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=8569738" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=gregmarr">
      <img class="avatar" src="https://avatars0.githubusercontent.com/u/8569738?s=40&amp;v=4" width="20" height="20" alt="@gregmarr" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=26778249" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=zerodefect">
      <img class="avatar" src="https://avatars2.githubusercontent.com/u/26778249?s=40&amp;v=4" width="20" height="20" alt="@zerodefect" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=5544610" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=Dobiasd">
      <img class="avatar" src="https://avatars0.githubusercontent.com/u/5544610?s=40&amp;v=4" width="20" height="20" alt="@Dobiasd" /> 
</a>    <a class="avatar-link" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=1706857" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/nlohmann/json/commits/develop/README.md?author=dan-42">
      <img class="avatar" src="https://avatars3.githubusercontent.com/u/1706857?s=40&amp;v=4" width="20" height="20" alt="@dan-42" /> 
</a>
    <button type="button" class="btn-link" data-toggle-for="blob_contributors_box">and others</button>

    </div>
  </div>



    <div class="file ">
      <div class="file-header">
  <div class="file-actions">


    <div class="BtnGroup">
      <a id="raw-url" class="btn btn-sm BtnGroup-item" href="/nlohmann/json/raw/develop/README.md">Raw</a>
        <a class="btn btn-sm js-update-url-with-hash BtnGroup-item" data-hotkey="b" href="/nlohmann/json/blame/develop/README.md">Blame</a>
      <a rel="nofollow" class="btn btn-sm BtnGroup-item" href="/nlohmann/json/commits/develop/README.md">History</a>
    </div>


        <button type="button" class="btn-octicon disabled tooltipped tooltipped-nw"
          aria-label="You must be signed in to make or propose changes">
          <svg class="octicon octicon-pencil" viewBox="0 0 14 16" version="1.1" width="14" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M0 12v3h3l8-8-3-3-8 8zm3 2H1v-2h1v1h1v1zm10.3-9.3L12 6 9 3l1.3-1.3a.996.996 0 0 1 1.41 0l1.59 1.59c.39.39.39 1.02 0 1.41z"/></svg>
        </button>
        <button type="button" class="btn-octicon btn-octicon-danger disabled tooltipped tooltipped-nw"
          aria-label="You must be signed in to make or propose changes">
          <svg class="octicon octicon-trashcan" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M11 2H9c0-.55-.45-1-1-1H5c-.55 0-1 .45-1 1H2c-.55 0-1 .45-1 1v1c0 .55.45 1 1 1v9c0 .55.45 1 1 1h7c.55 0 1-.45 1-1V5c.55 0 1-.45 1-1V3c0-.55-.45-1-1-1zm-1 12H3V5h1v8h1V5h1v8h1V5h1v8h1V5h1v9zm1-10H2V3h9v1z"/></svg>
        </button>
  </div>

  <div class="file-info">
      1268 lines (986 sloc)
      <span class="file-info-divider"></span>
    67.3 KB
  </div>
</div>

      
  <div id="readme" class="readme blob instapaper_body">
    <article class="markdown-body entry-content" itemprop="text"><p><a href="https://github.com/nlohmann/json/releases"><img src="https://raw.githubusercontent.com/nlohmann/json/master/doc/json.gif" alt="JSON for Modern C++" style="max-width:100%;"></a></p>
<p><a href="https://travis-ci.org/nlohmann/json" rel="nofollow"><img src="https://camo.githubusercontent.com/062bab570f5c0f60f166189c0e83160bf9c15da1/68747470733a2f2f7472617669732d63692e6f72672f6e6c6f686d616e6e2f6a736f6e2e7376673f6272616e63683d6d6173746572" alt="Build Status" data-canonical-src="https://travis-ci.org/nlohmann/json.svg?branch=master" style="max-width:100%;"></a>
<a href="https://ci.appveyor.com/project/nlohmann/json" rel="nofollow"><img src="https://camo.githubusercontent.com/0431ed1e58f0186dd67a419906cd5f33edd9fe41/68747470733a2f2f63692e6170707665796f722e636f6d2f6170692f70726f6a656374732f7374617475732f3161636233363678667967337179626b2f6272616e63682f646576656c6f703f7376673d74727565" alt="Build Status" data-canonical-src="https://ci.appveyor.com/api/projects/status/1acb366xfyg3qybk/branch/develop?svg=true" style="max-width:100%;"></a>
<a href="https://coveralls.io/r/nlohmann/json" rel="nofollow"><img src="https://camo.githubusercontent.com/167bd85bffbdf15a74a0ef79e6a342d1d1aabe71/68747470733a2f2f696d672e736869656c64732e696f2f636f766572616c6c732f6e6c6f686d616e6e2f6a736f6e2e737667" alt="Coverage Status" data-canonical-src="https://img.shields.io/coveralls/nlohmann/json.svg" style="max-width:100%;"></a>
<a href="https://scan.coverity.com/projects/nlohmann-json" rel="nofollow"><img src="https://camo.githubusercontent.com/bb0792bdccc05aa0c49a7cf4e8f275a270bff869/68747470733a2f2f7363616e2e636f7665726974792e636f6d2f70726f6a656374732f353535302f62616467652e737667" alt="Coverity Scan Build Status" data-canonical-src="https://scan.coverity.com/projects/5550/badge.svg" style="max-width:100%;"></a>
<a href="https://www.codacy.com/app/nlohmann/json?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=nlohmann/json&amp;utm_campaign=Badge_Grade" rel="nofollow"><img src="https://camo.githubusercontent.com/01c91aed99aad98d52a144e8483cbebf29f29d28/68747470733a2f2f6170692e636f646163792e636f6d2f70726f6a6563742f62616467652f47726164652f6633373332623333323765333433353861306539643166653966363631663038" alt="Codacy Badge" data-canonical-src="https://api.codacy.com/project/badge/Grade/f3732b3327e34358a0e9d1fe9f661f08" style="max-width:100%;"></a>
<a href="https://lgtm.com/projects/g/nlohmann/json/context:cpp" rel="nofollow"><img src="https://camo.githubusercontent.com/dc19be93a678f26c37f146800aa25b49c74dd447/68747470733a2f2f696d672e736869656c64732e696f2f6c67746d2f67726164652f6370702f672f6e6c6f686d616e6e2f6a736f6e2e7376673f6c6f676f3d6c67746d266c6f676f57696474683d3138" alt="Language grade: C/C++" data-canonical-src="https://img.shields.io/lgtm/grade/cpp/g/nlohmann/json.svg?logo=lgtm&amp;logoWidth=18" style="max-width:100%;"></a>
<a href="https://wandbox.org/permlink/TarF5pPn9NtHQjhf" rel="nofollow"><img src="https://camo.githubusercontent.com/57993e2729bc1a6482ef48e7c2fc3676a2df9c06/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f7472792d6f6e6c696e652d626c75652e737667" alt="Try online" data-canonical-src="https://img.shields.io/badge/try-online-blue.svg" style="max-width:100%;"></a>
<a href="http://nlohmann.github.io/json" rel="nofollow"><img src="https://camo.githubusercontent.com/ad89d109e98660d91ac965c4517e4f784132cdee/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f646f63732d646f787967656e2d626c75652e737667" alt="Documentation" data-canonical-src="https://img.shields.io/badge/docs-doxygen-blue.svg" style="max-width:100%;"></a>
<a href="https://raw.githubusercontent.com/nlohmann/json/master/LICENSE.MIT" rel="nofollow"><img src="https://camo.githubusercontent.com/890acbdcb87868b382af9a4b1fac507b9659d9bf/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f6c6963656e73652d4d49542d626c75652e737667" alt="GitHub license" data-canonical-src="https://img.shields.io/badge/license-MIT-blue.svg" style="max-width:100%;"></a>
<a href="https://github.com/nlohmann/json/releases"><img src="https://camo.githubusercontent.com/16fdc80c4217085da59920ac8898c82bfe907f13/68747470733a2f2f696d672e736869656c64732e696f2f6769746875622f72656c656173652f6e6c6f686d616e6e2f6a736f6e2e737667" alt="GitHub Releases" data-canonical-src="https://img.shields.io/github/release/nlohmann/json.svg" style="max-width:100%;"></a>
<a href="http://github.com/nlohmann/json/issues"><img src="https://camo.githubusercontent.com/7f829ff1c8b5c50d1717fd5363d7de35d11881c2/68747470733a2f2f696d672e736869656c64732e696f2f6769746875622f6973737565732f6e6c6f686d616e6e2f6a736f6e2e737667" alt="GitHub Issues" data-canonical-src="https://img.shields.io/github/issues/nlohmann/json.svg" style="max-width:100%;"></a>
<a href="http://isitmaintained.com/project/nlohmann/json" title="Average time to resolve an issue" rel="nofollow"><img src="https://camo.githubusercontent.com/8e6921366566ff9db7b685b83bfbfac17fe24091/687474703a2f2f697369746d61696e7461696e65642e636f6d2f62616467652f7265736f6c7574696f6e2f6e6c6f686d616e6e2f6a736f6e2e737667" alt="Average time to resolve an issue" data-canonical-src="http://isitmaintained.com/badge/resolution/nlohmann/json.svg" style="max-width:100%;"></a>
<a href="https://bestpractices.coreinfrastructure.org/projects/289" rel="nofollow"><img src="https://camo.githubusercontent.com/11a7141b8a8e15d5c44df845d08e30b2afb17e4a/68747470733a2f2f626573747072616374696365732e636f7265696e6672617374727563747572652e6f72672f70726f6a656374732f3238392f6261646765" alt="CII Best Practices" data-canonical-src="https://bestpractices.coreinfrastructure.org/projects/289/badge" style="max-width:100%;"></a></p>
<ul>
<li><a href="#design-goals">Design goals</a></li>
<li><a href="#integration">Integration</a>
<ul>
<li><a href="#cmake">CMake</a></li>
<li><a href="#package-managers">Package Managers</a></li>
</ul>
</li>
<li><a href="#examples">Examples</a>
<ul>
<li><a href="#json-as-first-class-data-type">JSON as first-class data type</a></li>
<li><a href="#serialization--deserialization">Serialization / Deserialization</a></li>
<li><a href="#stl-like-access">STL-like access</a></li>
<li><a href="#conversion-from-stl-containers">Conversion from STL containers</a></li>
<li><a href="#json-pointer-and-json-patch">JSON Pointer and JSON Patch</a></li>
<li><a href="#json-merge-patch">JSON Merge Patch</a></li>
<li><a href="#implicit-conversions">Implicit conversions</a></li>
<li><a href="#arbitrary-types-conversions">Conversions to/from arbitrary types</a></li>
<li><a href="#specializing-enum-conversion">Specializing enum conversion</a></li>
<li><a href="#binary-formats-bson-cbor-messagepack-and-ubjson">Binary formats (BSON, CBOR, MessagePack, and UBJSON)</a></li>
</ul>
</li>
<li><a href="#supported-compilers">Supported compilers</a></li>
<li><a href="#license">License</a></li>
<li><a href="#contact">Contact</a></li>
<li><a href="#thanks">Thanks</a></li>
<li><a href="#used-third-party-tools">Used third-party tools</a></li>
<li><a href="#projects-using-json-for-modern-c">Projects using JSON for Modern C++</a></li>
<li><a href="#notes">Notes</a></li>
<li><a href="#execute-unit-tests">Execute unit tests</a></li>
</ul>
<h2><a id="user-content-design-goals" class="anchor" aria-hidden="true" href="#design-goals"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Design goals</h2>
<p>There are myriads of <a href="http://json.org" rel="nofollow">JSON</a> libraries out there, and each may even have its reason to exist. Our class had these design goals:</p>
<ul>
<li>
<p><strong>Intuitive syntax</strong>. In languages such as Python, JSON feels like a first class data type. We used all the operator magic of modern C++ to achieve the same feeling in your code. Check out the <a href="#examples">examples below</a> and you'll know what I mean.</p>
</li>
<li>
<p><strong>Trivial integration</strong>. Our whole code consists of a single header file <a href="https://github.com/nlohmann/json/blob/develop/single_include/nlohmann/json.hpp"><code>json.hpp</code></a>. That's it. No library, no subproject, no dependencies, no complex build system. The class is written in vanilla C++11. All in all, everything should require no adjustment of your compiler flags or project settings.</p>
</li>
<li>
<p><strong>Serious testing</strong>. Our class is heavily <a href="https://github.com/nlohmann/json/tree/develop/test/src">unit-tested</a> and covers <a href="https://coveralls.io/r/nlohmann/json" rel="nofollow">100%</a> of the code, including all exceptional behavior. Furthermore, we checked with <a href="http://valgrind.org" rel="nofollow">Valgrind</a> and the <a href="https://clang.llvm.org/docs/index.html" rel="nofollow">Clang Sanitizers</a> that there are no memory leaks. <a href="https://github.com/google/oss-fuzz/tree/master/projects/json">Google OSS-Fuzz</a> additionally runs fuzz tests agains all parsers 24/7, effectively executing billions of tests so far. To maintain high quality, the project is following the <a href="https://bestpractices.coreinfrastructure.org/projects/289" rel="nofollow">Core Infrastructure Initiative (CII) best practices</a>.</p>
</li>
</ul>
<p>Other aspects were not so important to us:</p>
<ul>
<li>
<p><strong>Memory efficiency</strong>. Each JSON object has an overhead of one pointer (the maximal size of a union) and one enumeration element (1 byte). The default generalization uses the following C++ data types: <code>std::string</code> for strings, <code>int64_t</code>, <code>uint64_t</code> or <code>double</code> for numbers, <code>std::map</code> for objects, <code>std::vector</code> for arrays, and <code>bool</code> for Booleans. However, you can template the generalized class <code>basic_json</code> to your needs.</p>
</li>
<li>
<p><strong>Speed</strong>. There are certainly <a href="https://github.com/miloyip/nativejson-benchmark#parsing-time">faster JSON libraries</a> out there. However, if your goal is to speed up your development by adding JSON support with a single header, then this library is the way to go. If you know how to use a <code>std::vector</code> or <code>std::map</code>, you are already set.</p>
</li>
</ul>
<p>See the <a href="https://github.com/nlohmann/json/blob/master/.github/CONTRIBUTING.md#please-dont">contribution guidelines</a> for more information.</p>
<h2><a id="user-content-integration" class="anchor" aria-hidden="true" href="#integration"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Integration</h2>
<p><a href="https://github.com/nlohmann/json/blob/develop/single_include/nlohmann/json.hpp"><code>json.hpp</code></a> is the single required file in <code>single_include/nlohmann</code> or <a href="https://github.com/nlohmann/json/releases">released here</a>. You need to add</p>
<div class="highlight highlight-source-c++"><pre>#<span class="pl-k">include</span> <span class="pl-s"><span class="pl-pds">&lt;</span>nlohmann/json.hpp<span class="pl-pds">&gt;</span></span>

<span class="pl-c"><span class="pl-c">//</span> for convenience</span>
<span class="pl-k">using</span> json = nlohmann::json;</pre></div>
<p>to the files you want to process JSON and set the necessary switches to enable C++11 (e.g., <code>-std=c++11</code> for GCC and Clang).</p>
<p>You can further use file <a href="https://github.com/nlohmann/json/blob/develop/include/nlohmann/json_fwd.hpp"><code>include/nlohmann/json_fwd.hpp</code></a> for forward-declarations. The installation of json_fwd.hpp (as part of cmake's install step), can be achieved by setting <code>-DJSON_MultipleHeaders=ON</code>.</p>
<h3><a id="user-content-cmake" class="anchor" aria-hidden="true" href="#cmake"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>CMake</h3>
<p>You can also use the <code>nlohmann_json::nlohmann_json</code> interface target in CMake.  This target populates the appropriate usage requirements for <code>INTERFACE_INCLUDE_DIRECTORIES</code> to point to the appropriate include directories and <code>INTERFACE_COMPILE_FEATURES</code> for the necessary C++11 flags.</p>
<h4><a id="user-content-external" class="anchor" aria-hidden="true" href="#external"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>External</h4>
<p>To use this library from a CMake project, you can locate it directly with <code>find_package()</code> and use the namespaced imported target from the generated package configuration:</p>
<div class="highlight highlight-source-cmake"><pre><span class="pl-c"><span class="pl-c">#</span> CMakeLists.txt</span>
<span class="pl-c1">find_package</span>(nlohmann_json 3.2.0 <span class="pl-k">REQUIRED</span>)
...
<span class="pl-c1">add_library</span>(foo ...)
...
<span class="pl-c1">target_link_libraries</span>(foo <span class="pl-k">PRIVATE</span> nlohmann_json::nlohmann_json)</pre></div>
<p>The package configuration file, <code>nlohmann_jsonConfig.cmake</code>, can be used either from an install tree or directly out of the build tree.</p>
<h4><a id="user-content-embedded" class="anchor" aria-hidden="true" href="#embedded"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Embedded</h4>
<p>To embed the library directly into an existing CMake project, place the entire source tree in a subdirectory and call <code>add_subdirectory()</code> in your <code>CMakeLists.txt</code> file:</p>
<div class="highlight highlight-source-cmake"><pre><span class="pl-c"><span class="pl-c">#</span> Typically you don't care so much for a third party library's tests to be</span>
<span class="pl-c"><span class="pl-c">#</span> run from your own project's code.</span>
<span class="pl-c1">set</span>(JSON_BuildTests <span class="pl-k">OFF</span> <span class="pl-k">CACHE</span> INTERNAL <span class="pl-s">""</span>)

<span class="pl-c"><span class="pl-c">#</span> Don't use include(nlohmann_json/CMakeLists.txt) since that carries with it</span>
<span class="pl-c"><span class="pl-c">#</span> inintended consequences that will break the build.  It's generally</span>
<span class="pl-c"><span class="pl-c">#</span> discouraged (although not necessarily well documented as such) to use</span>
<span class="pl-c"><span class="pl-c">#</span> include(...) for pulling in other CMake projects anyways.</span>
<span class="pl-c1">add_subdirectory</span>(nlohmann_json)
...
<span class="pl-c1">add_library</span>(foo ...)
...
<span class="pl-c1">target_link_libraries</span>(foo <span class="pl-k">PRIVATE</span> nlohmann_json::nlohmann_json)</pre></div>
<h4><a id="user-content-supporting-both" class="anchor" aria-hidden="true" href="#supporting-both"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Supporting Both</h4>
<p>To allow your project to support either an externally supplied or an embedded JSON library, you can use a pattern akin to the following:</p>
<div class="highlight highlight-source-cmake"><pre><span class="pl-c"><span class="pl-c">#</span> Top level CMakeLists.txt</span>
<span class="pl-c1">project</span>(FOO)
...
<span class="pl-c1">option</span>(FOO_USE_EXTERNAL_JSON <span class="pl-s">"Use an external JSON library"</span> <span class="pl-k">OFF</span>)
...
<span class="pl-c1">add_subdirectory</span>(thirdparty)
...
<span class="pl-c1">add_library</span>(foo ...)
...
<span class="pl-c"><span class="pl-c">#</span> Note that the namespaced target will always be available regardless of the</span>
<span class="pl-c"><span class="pl-c">#</span> import method</span>
<span class="pl-c1">target_link_libraries</span>(foo <span class="pl-k">PRIVATE</span> nlohmann_json::nlohmann_json)</pre></div>
<div class="highlight highlight-source-cmake"><pre><span class="pl-c"><span class="pl-c">#</span> thirdparty/CMakeLists.txt</span>
...
<span class="pl-k">if</span>(FOO_USE_EXTERNAL_JSON)
  <span class="pl-c1">find_package</span>(nlohmann_json 3.2.0 <span class="pl-k">REQUIRED</span>)
<span class="pl-k">else</span>()
  <span class="pl-c1">set</span>(JSON_BuildTests <span class="pl-k">OFF</span> <span class="pl-k">CACHE</span> INTERNAL <span class="pl-s">""</span>)
  <span class="pl-c1">add_subdirectory</span>(nlohmann_json)
<span class="pl-k">endif</span>()
...</pre></div>
<p><code>thirdparty/nlohmann_json</code> is then a complete copy of this source tree.</p>
<h3><a id="user-content-package-managers" class="anchor" aria-hidden="true" href="#package-managers"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Package Managers</h3>
<p><g-emoji class="g-emoji" alias="beer" fallback-src="https://assets-cdn.github.com/images/icons/emoji/unicode/1f37a.png">🍺</g-emoji> If you are using OS X and <a href="http://brew.sh" rel="nofollow">Homebrew</a>, just type <code>brew tap nlohmann/json</code> and <code>brew install nlohmann_json</code> and you're set. If you want the bleeding edge rather than the latest release, use <code>brew install nlohmann_json --HEAD</code>.</p>
<p>If you are using the <a href="http://mesonbuild.com" rel="nofollow">Meson Build System</a>, then you can get a wrap file by downloading it from <a href="https://wrapdb.mesonbuild.com/nlohmann_json" rel="nofollow">Meson WrapDB</a>, or simply use <code>meson wrap install nlohmann_json</code>.</p>
<p>If you are using <a href="https://www.conan.io/" rel="nofollow">Conan</a> to manage your dependencies, merely add <code>jsonformoderncpp/x.y.z@vthiery/stable</code> to your <code>conanfile.py</code>'s requires, where <code>x.y.z</code> is the release version you want to use. Please file issues <a href="https://github.com/vthiery/conan-jsonformoderncpp/issues">here</a> if you experience problems with the packages.</p>
<p>If you are using <a href="https://www.spack.io/" rel="nofollow">Spack</a> to manage your dependencies, you can use the <code>nlohmann_json</code> package. Please see the <a href="https://github.com/spack/spack">spack project</a> for any issues regarding the packaging.</p>
<p>If you are using <a href="https://github.com/ruslo/hunter/">hunter</a> on your project for external dependencies, then you can use the <a href="https://docs.hunter.sh/en/latest/packages/pkg/nlohmann_json.html" rel="nofollow">nlohmann_json package</a>. Please see the hunter project for any issues regarding the packaging.</p>
<p>If you are using <a href="https://buckaroo.pm" rel="nofollow">Buckaroo</a>, you can install this library's module with <code>buckaroo install nlohmann/json</code>. Please file issues <a href="https://github.com/LoopPerfect/buckaroo-recipes/issues/new?title=nlohmann/nlohmann/json">here</a>.</p>
<p>If you are using <a href="https://github.com/Microsoft/vcpkg/">vcpkg</a> on your project for external dependencies, then you can use the <a href="https://github.com/Microsoft/vcpkg/tree/master/ports/nlohmann-json">nlohmann-json package</a>. Please see the vcpkg project for any issues regarding the packaging.</p>
<p>If you are using <a href="http://cget.readthedocs.io/en/latest/" rel="nofollow">cget</a>, you can install the latest development version with <code>cget install nlohmann/json</code>. A specific version can be installed with <code>cget install nlohmann/json@v3.1.0</code>. Also, the multiple header version can be installed by adding the <code>-DJSON_MultipleHeaders=ON</code> flag (i.e., <code>cget install nlohmann/json -DJSON_MultipleHeaders=ON</code>).</p>
<p>If you are using <a href="https://cocoapods.org" rel="nofollow">CocoaPods</a>, you can use the library by adding pod <code>"nlohmann_json", '~&gt;3.1.2'</code> to your podfile (see <a href="https://bitbucket.org/benman/nlohmann_json-cocoapod/src/master/" rel="nofollow">an example</a>). Please file issues <a href="https://bitbucket.org/benman/nlohmann_json-cocoapod/issues?status=new&amp;status=open" rel="nofollow">here</a>.</p>
<h2><a id="user-content-examples" class="anchor" aria-hidden="true" href="#examples"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Examples</h2>
<p>Beside the examples below, you may want to check the <a href="https://nlohmann.github.io/json/" rel="nofollow">documentation</a> where each function contains a separate code example (e.g., check out <a href="https://nlohmann.github.io/json/classnlohmann_1_1basic__json_a5338e282d1d02bed389d852dd670d98d.html#a5338e282d1d02bed389d852dd670d98d" rel="nofollow"><code>emplace()</code></a>). All <a href="https://github.com/nlohmann/json/tree/develop/doc/examples">example files</a> can be compiled and executed on their own (e.g., file <a href="https://github.com/nlohmann/json/blob/develop/doc/examples/emplace.cpp">emplace.cpp</a>).</p>
<h3><a id="user-content-json-as-first-class-data-type" class="anchor" aria-hidden="true" href="#json-as-first-class-data-type"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>JSON as first-class data type</h3>
<p>Here are some examples to give you an idea how to use the class.</p>
<p>Assume you want to create the JSON object</p>
<div class="highlight highlight-source-json"><pre>{
  <span class="pl-s"><span class="pl-pds">"</span>pi<span class="pl-pds">"</span></span>: <span class="pl-c1">3.141</span>,
  <span class="pl-s"><span class="pl-pds">"</span>happy<span class="pl-pds">"</span></span>: <span class="pl-c1">true</span>,
  <span class="pl-s"><span class="pl-pds">"</span>name<span class="pl-pds">"</span></span>: <span class="pl-s"><span class="pl-pds">"</span>Niels<span class="pl-pds">"</span></span>,
  <span class="pl-s"><span class="pl-pds">"</span>nothing<span class="pl-pds">"</span></span>: <span class="pl-c1">null</span>,
  <span class="pl-s"><span class="pl-pds">"</span>answer<span class="pl-pds">"</span></span>: {
    <span class="pl-s"><span class="pl-pds">"</span>everything<span class="pl-pds">"</span></span>: <span class="pl-c1">42</span>
  },
  <span class="pl-s"><span class="pl-pds">"</span>list<span class="pl-pds">"</span></span>: [<span class="pl-c1">1</span>, <span class="pl-c1">0</span>, <span class="pl-c1">2</span>],
  <span class="pl-s"><span class="pl-pds">"</span>object<span class="pl-pds">"</span></span>: {
    <span class="pl-s"><span class="pl-pds">"</span>currency<span class="pl-pds">"</span></span>: <span class="pl-s"><span class="pl-pds">"</span>USD<span class="pl-pds">"</span></span>,
    <span class="pl-s"><span class="pl-pds">"</span>value<span class="pl-pds">"</span></span>: <span class="pl-c1">42.99</span>
  }
}</pre></div>
<p>With this library, you could write:</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> create an empty structure (null)</span>
json j;

<span class="pl-c"><span class="pl-c">//</span> add a number that is stored as double (note the implicit conversion of j to an object)</span>
j[<span class="pl-s"><span class="pl-pds">"</span>pi<span class="pl-pds">"</span></span>] = <span class="pl-c1">3.141</span>;

<span class="pl-c"><span class="pl-c">//</span> add a Boolean that is stored as bool</span>
j[<span class="pl-s"><span class="pl-pds">"</span>happy<span class="pl-pds">"</span></span>] = <span class="pl-c1">true</span>;

<span class="pl-c"><span class="pl-c">//</span> add a string that is stored as std::string</span>
j[<span class="pl-s"><span class="pl-pds">"</span>name<span class="pl-pds">"</span></span>] = <span class="pl-s"><span class="pl-pds">"</span>Niels<span class="pl-pds">"</span></span>;

<span class="pl-c"><span class="pl-c">//</span> add another null object by passing nullptr</span>
j[<span class="pl-s"><span class="pl-pds">"</span>nothing<span class="pl-pds">"</span></span>] = <span class="pl-c1">nullptr</span>;

<span class="pl-c"><span class="pl-c">//</span> add an object inside the object</span>
j[<span class="pl-s"><span class="pl-pds">"</span>answer<span class="pl-pds">"</span></span>][<span class="pl-s"><span class="pl-pds">"</span>everything<span class="pl-pds">"</span></span>] = <span class="pl-c1">42</span>;

<span class="pl-c"><span class="pl-c">//</span> add an array that is stored as std::vector (using an initializer list)</span>
j[<span class="pl-s"><span class="pl-pds">"</span>list<span class="pl-pds">"</span></span>] = { <span class="pl-c1">1</span>, <span class="pl-c1">0</span>, <span class="pl-c1">2</span> };

<span class="pl-c"><span class="pl-c">//</span> add another object (using an initializer list of pairs)</span>
j[<span class="pl-s"><span class="pl-pds">"</span>object<span class="pl-pds">"</span></span>] = { {<span class="pl-s"><span class="pl-pds">"</span>currency<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>USD<span class="pl-pds">"</span></span>}, {<span class="pl-s"><span class="pl-pds">"</span>value<span class="pl-pds">"</span></span>, <span class="pl-c1">42.99</span>} };

<span class="pl-c"><span class="pl-c">//</span> instead, you could also write (which looks very similar to the JSON above)</span>
json j2 = {
  {<span class="pl-s"><span class="pl-pds">"</span>pi<span class="pl-pds">"</span></span>, <span class="pl-c1">3.141</span>},
  {<span class="pl-s"><span class="pl-pds">"</span>happy<span class="pl-pds">"</span></span>, <span class="pl-c1">true</span>},
  {<span class="pl-s"><span class="pl-pds">"</span>name<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>Niels<span class="pl-pds">"</span></span>},
  {<span class="pl-s"><span class="pl-pds">"</span>nothing<span class="pl-pds">"</span></span>, <span class="pl-c1">nullptr</span>},
  {<span class="pl-s"><span class="pl-pds">"</span>answer<span class="pl-pds">"</span></span>, {
    {<span class="pl-s"><span class="pl-pds">"</span>everything<span class="pl-pds">"</span></span>, <span class="pl-c1">42</span>}
  }},
  {<span class="pl-s"><span class="pl-pds">"</span>list<span class="pl-pds">"</span></span>, {<span class="pl-c1">1</span>, <span class="pl-c1">0</span>, <span class="pl-c1">2</span>}},
  {<span class="pl-s"><span class="pl-pds">"</span>object<span class="pl-pds">"</span></span>, {
    {<span class="pl-s"><span class="pl-pds">"</span>currency<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>USD<span class="pl-pds">"</span></span>},
    {<span class="pl-s"><span class="pl-pds">"</span>value<span class="pl-pds">"</span></span>, <span class="pl-c1">42.99</span>}
  }}
};</pre></div>
<p>Note that in all these cases, you never need to "tell" the compiler which JSON value type you want to use. If you want to be explicit or express some edge cases, the functions <a href="https://nlohmann.github.io/json/classnlohmann_1_1basic__json_aa80485befaffcadaa39965494e0b4d2e.html#aa80485befaffcadaa39965494e0b4d2e" rel="nofollow"><code>json::array</code></a> and <a href="https://nlohmann.github.io/json/classnlohmann_1_1basic__json_aa13f7c0615867542ce80337cbcf13ada.html#aa13f7c0615867542ce80337cbcf13ada" rel="nofollow"><code>json::object</code></a> will help:</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> a way to express the empty array []</span>
json empty_array_explicit = json::array();

<span class="pl-c"><span class="pl-c">//</span> ways to express the empty object {}</span>
json empty_object_implicit = json({});
json empty_object_explicit = json::object();

<span class="pl-c"><span class="pl-c">//</span> a way to express an _array_ of key/value pairs [["currency", "USD"], ["value", 42.99]]</span>
json array_not_object = json::array({ {<span class="pl-s"><span class="pl-pds">"</span>currency<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>USD<span class="pl-pds">"</span></span>}, {<span class="pl-s"><span class="pl-pds">"</span>value<span class="pl-pds">"</span></span>, <span class="pl-c1">42.99</span>} });</pre></div>
<h3><a id="user-content-serialization--deserialization" class="anchor" aria-hidden="true" href="#serialization--deserialization"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Serialization / Deserialization</h3>
<h4><a id="user-content-tofrom-strings" class="anchor" aria-hidden="true" href="#tofrom-strings"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>To/from strings</h4>
<p>You can create a JSON value (deserialization) by appending <code>_json</code> to a string literal:</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> create object from string literal</span>
json j = <span class="pl-s"><span class="pl-pds">"</span>{ <span class="pl-cce">\"</span>happy<span class="pl-cce">\"</span>: true, <span class="pl-cce">\"</span>pi<span class="pl-cce">\"</span>: 3.141 }<span class="pl-pds">"</span></span>_json;

<span class="pl-c"><span class="pl-c">//</span> or even nicer with a raw string literal</span>
<span class="pl-k">auto</span> j2 = <span class="pl-s"><span class="pl-pds">R"(</span></span>
<span class="pl-s">  {</span>
<span class="pl-s">    "happy": true,</span>
<span class="pl-s">    "pi": 3.141</span>
<span class="pl-s">  }</span>
<span class="pl-s"><span class="pl-pds">)"</span></span>_json;</pre></div>
<p>Note that without appending the <code>_json</code> suffix, the passed string literal is not parsed, but just used as JSON string value. That is, <code>json j = "{ \"happy\": true, \"pi\": 3.141 }"</code> would just store the string <code>"{ "happy": true, "pi": 3.141 }"</code> rather than parsing the actual object.</p>
<p>The above example can also be expressed explicitly using <a href="https://nlohmann.github.io/json/classnlohmann_1_1basic__json_aa9676414f2e36383c4b181fe856aa3c0.html#aa9676414f2e36383c4b181fe856aa3c0" rel="nofollow"><code>json::parse()</code></a>:</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> parse explicitly</span>
<span class="pl-k">auto</span> j3 = json::parse(<span class="pl-s"><span class="pl-pds">"</span>{ <span class="pl-cce">\"</span>happy<span class="pl-cce">\"</span>: true, <span class="pl-cce">\"</span>pi<span class="pl-cce">\"</span>: 3.141 }<span class="pl-pds">"</span></span>);</pre></div>
<p>You can also get a string representation of a JSON value (serialize):</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> explicit conversion to string</span>
std::string s = j.dump();    <span class="pl-c"><span class="pl-c">//</span> {\"happy\":true,\"pi\":3.141}</span>

<span class="pl-c"><span class="pl-c">//</span> serialization with pretty printing</span>
<span class="pl-c"><span class="pl-c">//</span> pass in the amount of spaces to indent</span>
std::cout &lt;&lt; j.dump(<span class="pl-c1">4</span>) &lt;&lt; std::endl;
<span class="pl-c"><span class="pl-c">//</span> {</span>
<span class="pl-c"><span class="pl-c">//</span>     "happy": true,</span>
<span class="pl-c"><span class="pl-c">//</span>     "pi": 3.141</span>
<span class="pl-c"><span class="pl-c">//</span> }</span></pre></div>
<p>Note the difference between serialization and assignment:</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> store a string in a JSON value</span>
json j_string = <span class="pl-s"><span class="pl-pds">"</span>this is a string<span class="pl-pds">"</span></span>;

<span class="pl-c"><span class="pl-c">//</span> retrieve the string value (implicit JSON to std::string conversion)</span>
std::string cpp_string = j_string;
<span class="pl-c"><span class="pl-c">//</span> retrieve the string value (explicit JSON to std::string conversion)</span>
<span class="pl-k">auto</span> cpp_string2 = j_string.get&lt;std::string&gt;();
<span class="pl-c"><span class="pl-c">//</span> retrieve the string value (alternative explicit JSON to std::string conversion)</span>
std::string cpp_string3;
j_string.get_to(cpp_string3);

<span class="pl-c"><span class="pl-c">//</span> retrieve the serialized value (explicit JSON serialization)</span>
std::string serialized_string = j_string.dump();

<span class="pl-c"><span class="pl-c">//</span> output of original string</span>
std::cout &lt;&lt; cpp_string &lt;&lt; <span class="pl-s"><span class="pl-pds">"</span> == <span class="pl-pds">"</span></span> &lt;&lt; cpp_string2 &lt;&lt; <span class="pl-s"><span class="pl-pds">"</span> == <span class="pl-pds">"</span></span> &lt;&lt; cpp_string3 &lt;&lt; <span class="pl-s"><span class="pl-pds">"</span> == <span class="pl-pds">"</span></span> &lt;&lt; j_string.get&lt;std::string&gt;() &lt;&lt; <span class="pl-s"><span class="pl-pds">'</span><span class="pl-cce">\n</span><span class="pl-pds">'</span></span>;
<span class="pl-c"><span class="pl-c">//</span> output of serialized value</span>
std::cout &lt;&lt; j_string &lt;&lt; <span class="pl-s"><span class="pl-pds">"</span> == <span class="pl-pds">"</span></span> &lt;&lt; serialized_string &lt;&lt; std::endl;</pre></div>
<p><a href="https://nlohmann.github.io/json/classnlohmann_1_1basic__json_a50ec80b02d0f3f51130d4abb5d1cfdc5.html#a50ec80b02d0f3f51130d4abb5d1cfdc5" rel="nofollow"><code>.dump()</code></a> always returns the serialized value, and <a href="https://nlohmann.github.io/json/classnlohmann_1_1basic__json_a16f9445f7629f634221a42b967cdcd43.html#a16f9445f7629f634221a42b967cdcd43" rel="nofollow"><code>.get&lt;std::string&gt;()</code></a> returns the originally stored string value.</p>
<p>Note the library only supports UTF-8. When you store strings with different encodings in the library, calling <a href="https://nlohmann.github.io/json/classnlohmann_1_1basic__json_a50ec80b02d0f3f51130d4abb5d1cfdc5.html#a50ec80b02d0f3f51130d4abb5d1cfdc5" rel="nofollow"><code>dump()</code></a> may throw an exception unless <code>json::error_handler_t::replace</code> or <code>json::error_handler_t::ignore</code> are used as error handlers.</p>
<h4><a id="user-content-tofrom-streams-eg-files-string-streams" class="anchor" aria-hidden="true" href="#tofrom-streams-eg-files-string-streams"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>To/from streams (e.g. files, string streams)</h4>
<p>You can also use streams to serialize and deserialize:</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> deserialize from standard input</span>
json j;
std::cin &gt;&gt; j;

<span class="pl-c"><span class="pl-c">//</span> serialize to standard output</span>
std::cout &lt;&lt; j;

<span class="pl-c"><span class="pl-c">//</span> the setw manipulator was overloaded to set the indentation for pretty printing</span>
std::cout &lt;&lt; std::setw(<span class="pl-c1">4</span>) &lt;&lt; j &lt;&lt; std::endl;</pre></div>
<p>These operators work for any subclasses of <code>std::istream</code> or <code>std::ostream</code>. Here is the same example with files:</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> read a JSON file</span>
std::ifstream <span class="pl-en">i</span>(<span class="pl-s"><span class="pl-pds">"</span>file.json<span class="pl-pds">"</span></span>);
json j;
i &gt;&gt; j;

<span class="pl-c"><span class="pl-c">//</span> write prettified JSON to another file</span>
std::ofstream <span class="pl-en">o</span>(<span class="pl-s"><span class="pl-pds">"</span>pretty.json<span class="pl-pds">"</span></span>);
o &lt;&lt; std::setw(<span class="pl-c1">4</span>) &lt;&lt; j &lt;&lt; std::endl;</pre></div>
<p>Please note that setting the exception bit for <code>failbit</code> is inappropriate for this use case. It will result in program termination due to the <code>noexcept</code> specifier in use.</p>
<h4><a id="user-content-read-from-iterator-range" class="anchor" aria-hidden="true" href="#read-from-iterator-range"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Read from iterator range</h4>
<p>You can also parse JSON from an iterator range; that is, from any container accessible by iterators whose content is stored as contiguous byte sequence, for instance a <code>std::vector&lt;std::uint8_t&gt;</code>:</p>
<div class="highlight highlight-source-c++"><pre>std::vector&lt;std::<span class="pl-c1">uint8_t</span>&gt; v = {<span class="pl-s"><span class="pl-pds">'</span>t<span class="pl-pds">'</span></span>, <span class="pl-s"><span class="pl-pds">'</span>r<span class="pl-pds">'</span></span>, <span class="pl-s"><span class="pl-pds">'</span>u<span class="pl-pds">'</span></span>, <span class="pl-s"><span class="pl-pds">'</span>e<span class="pl-pds">'</span></span>};
json j = json::parse(v.begin(), v.end());</pre></div>
<p>You may leave the iterators for the range [begin, end):</p>
<div class="highlight highlight-source-c++"><pre>std::vector&lt;std::<span class="pl-c1">uint8_t</span>&gt; v = {<span class="pl-s"><span class="pl-pds">'</span>t<span class="pl-pds">'</span></span>, <span class="pl-s"><span class="pl-pds">'</span>r<span class="pl-pds">'</span></span>, <span class="pl-s"><span class="pl-pds">'</span>u<span class="pl-pds">'</span></span>, <span class="pl-s"><span class="pl-pds">'</span>e<span class="pl-pds">'</span></span>};
json j = json::parse(v);</pre></div>
<h4><a id="user-content-sax-interface" class="anchor" aria-hidden="true" href="#sax-interface"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>SAX interface</h4>
<p>The library uses a SAX-like interface with the following functions:</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> called when null is parsed</span>
<span class="pl-k">bool</span> <span class="pl-en">null</span>();

<span class="pl-c"><span class="pl-c">//</span> called when a boolean is parsed; value is passed</span>
<span class="pl-k">bool</span> <span class="pl-en">boolean</span>(<span class="pl-k">bool</span> val);

<span class="pl-c"><span class="pl-c">//</span> called when a signed or unsigned integer number is parsed; value is passed</span>
<span class="pl-k">bool</span> <span class="pl-en">number_integer</span>(<span class="pl-c1">number_integer_t</span> val);
<span class="pl-k">bool</span> <span class="pl-en">number_unsigned</span>(<span class="pl-c1">number_unsigned_t</span> val);

<span class="pl-c"><span class="pl-c">//</span> called when a floating-point number is parsed; value and original string is passed</span>
<span class="pl-k">bool</span> <span class="pl-en">number_float</span>(<span class="pl-c1">number_float_t</span> val, <span class="pl-k">const</span> <span class="pl-c1">string_t</span>&amp; s);

<span class="pl-c"><span class="pl-c">//</span> called when a string is parsed; value is passed and can be safely moved away</span>
<span class="pl-k">bool</span> <span class="pl-en">string</span>(<span class="pl-c1">string_t</span>&amp; val);

<span class="pl-c"><span class="pl-c">//</span> called when an object or array begins or ends, resp. The number of elements is passed (or -1 if not known)</span>
<span class="pl-k">bool</span> <span class="pl-en">start_object</span>(std::<span class="pl-c1">size_t</span> elements);
<span class="pl-k">bool</span> <span class="pl-en">end_object</span>();
<span class="pl-k">bool</span> <span class="pl-en">start_array</span>(std::<span class="pl-c1">size_t</span> elements);
<span class="pl-k">bool</span> <span class="pl-en">end_array</span>();
<span class="pl-c"><span class="pl-c">//</span> called when an object key is parsed; value is passed and can be safely moved away</span>
<span class="pl-k">bool</span> <span class="pl-en">key</span>(<span class="pl-c1">string_t</span>&amp; val);

<span class="pl-c"><span class="pl-c">//</span> called when a parse error occurs; byte position, the last token, and an exception is passed</span>
<span class="pl-k">bool</span> <span class="pl-en">parse_error</span>(std::<span class="pl-c1">size_t</span> position, <span class="pl-k">const</span> std::string&amp; last_token, <span class="pl-k">const</span> detail::exception&amp; ex);</pre></div>
<p>The return value of each function determines whether parsing should proceed.</p>
<p>To implement your own SAX handler, proceed as follows:</p>
<ol>
<li>Implement the SAX interface in a class. You can use class <code>nlohmann::json_sax&lt;json&gt;</code> as base class, but you can also use any class where the functions described above are implemented and public.</li>
<li>Create an object of your SAX interface class, e.g. <code>my_sax</code>.</li>
<li>Call <code>bool json::sax_parse(input, &amp;my_sax)</code>; where the first parameter can be any input like a string or an input stream and the second parameter is a pointer to your SAX interface.</li>
</ol>
<p>Note the <code>sax_parse</code> function only returns a <code>bool</code> indicating the result of the last executed SAX event. It does not return a  <code>json</code> value - it is up to you to decide what to do with the SAX events. Furthermore, no exceptions are thrown in case of a parse error - it is up to you what to do with the exception object passed to your <code>parse_error</code> implementation. Internally, the SAX interface is used for the DOM parser (class <code>json_sax_dom_parser</code>) as well as the acceptor (<code>json_sax_acceptor</code>), see file <a href="https://github.com/nlohmann/json/blob/develop/include/nlohmann/detail/input/json_sax.hpp"><code>json_sax.hpp</code></a>.</p>
<h3><a id="user-content-stl-like-access" class="anchor" aria-hidden="true" href="#stl-like-access"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>STL-like access</h3>
<p>We designed the JSON class to behave just like an STL container. In fact, it satisfies the <a href="https://en.cppreference.com/w/cpp/named_req/ReversibleContainer" rel="nofollow"><strong>ReversibleContainer</strong></a> requirement.</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> create an array using push_back</span>
json j;
j.push_back(<span class="pl-s"><span class="pl-pds">"</span>foo<span class="pl-pds">"</span></span>);
j.push_back(<span class="pl-c1">1</span>);
j.push_back(<span class="pl-c1">true</span>);

<span class="pl-c"><span class="pl-c">//</span> also use emplace_back</span>
j.emplace_back(<span class="pl-c1">1.78</span>);

<span class="pl-c"><span class="pl-c">//</span> iterate the array</span>
<span class="pl-k">for</span> (json::iterator it = j.begin(); it != j.end(); ++it) {
  std::cout &lt;&lt; *it &lt;&lt; <span class="pl-s"><span class="pl-pds">'</span><span class="pl-cce">\n</span><span class="pl-pds">'</span></span>;
}

<span class="pl-c"><span class="pl-c">//</span> range-based for</span>
<span class="pl-k">for</span> (<span class="pl-k">auto</span>&amp; element : j) {
  std::cout &lt;&lt; element &lt;&lt; <span class="pl-s"><span class="pl-pds">'</span><span class="pl-cce">\n</span><span class="pl-pds">'</span></span>;
}

<span class="pl-c"><span class="pl-c">//</span> getter/setter</span>
<span class="pl-k">const</span> std::string tmp = j[<span class="pl-c1">0</span>];
j[<span class="pl-c1">1</span>] = <span class="pl-c1">42</span>;
<span class="pl-k">bool</span> foo = j.at(<span class="pl-c1">2</span>);

<span class="pl-c"><span class="pl-c">//</span> comparison</span>
j == <span class="pl-s"><span class="pl-pds">"</span>[<span class="pl-cce">\"</span>foo<span class="pl-cce">\"</span>, 1, true]<span class="pl-pds">"</span></span>_json;  <span class="pl-c"><span class="pl-c">//</span> true</span>

<span class="pl-c"><span class="pl-c">//</span> other stuff</span>
j.size();     <span class="pl-c"><span class="pl-c">//</span> 3 entries</span>
j.empty();    <span class="pl-c"><span class="pl-c">//</span> false</span>
j.type();     <span class="pl-c"><span class="pl-c">//</span> json::value_t::array</span>
j.clear();    <span class="pl-c"><span class="pl-c">//</span> the array is empty again</span>

<span class="pl-c"><span class="pl-c">//</span> convenience type checkers</span>
j.is_null();
j.is_boolean();
j.is_number();
j.is_object();
j.is_array();
j.is_string();

<span class="pl-c"><span class="pl-c">//</span> create an object</span>
json o;
o[<span class="pl-s"><span class="pl-pds">"</span>foo<span class="pl-pds">"</span></span>] = <span class="pl-c1">23</span>;
o[<span class="pl-s"><span class="pl-pds">"</span>bar<span class="pl-pds">"</span></span>] = <span class="pl-c1">false</span>;
o[<span class="pl-s"><span class="pl-pds">"</span>baz<span class="pl-pds">"</span></span>] = <span class="pl-c1">3.141</span>;

<span class="pl-c"><span class="pl-c">//</span> also use emplace</span>
o.emplace(<span class="pl-s"><span class="pl-pds">"</span>weather<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>sunny<span class="pl-pds">"</span></span>);

<span class="pl-c"><span class="pl-c">//</span> special iterator member functions for objects</span>
<span class="pl-k">for</span> (json::iterator it = o.begin(); it != o.end(); ++it) {
  std::cout &lt;&lt; it.<span class="pl-c1">key</span>() &lt;&lt; <span class="pl-s"><span class="pl-pds">"</span> : <span class="pl-pds">"</span></span> &lt;&lt; it.<span class="pl-c1">value</span>() &lt;&lt; <span class="pl-s"><span class="pl-pds">"</span><span class="pl-cce">\n</span><span class="pl-pds">"</span></span>;
}

<span class="pl-c"><span class="pl-c">//</span> find an entry</span>
<span class="pl-k">if</span> (o.find(<span class="pl-s"><span class="pl-pds">"</span>foo<span class="pl-pds">"</span></span>) != o.end()) {
  <span class="pl-c"><span class="pl-c">//</span> there is an entry with key "foo"</span>
}

<span class="pl-c"><span class="pl-c">//</span> or simpler using count()</span>
<span class="pl-k">int</span> foo_present = o.count(<span class="pl-s"><span class="pl-pds">"</span>foo<span class="pl-pds">"</span></span>); <span class="pl-c"><span class="pl-c">//</span> 1</span>
<span class="pl-k">int</span> fob_present = o.count(<span class="pl-s"><span class="pl-pds">"</span>fob<span class="pl-pds">"</span></span>); <span class="pl-c"><span class="pl-c">//</span> 0</span>

<span class="pl-c"><span class="pl-c">//</span> delete an entry</span>
o.erase(<span class="pl-s"><span class="pl-pds">"</span>foo<span class="pl-pds">"</span></span>);</pre></div>
<h3><a id="user-content-conversion-from-stl-containers" class="anchor" aria-hidden="true" href="#conversion-from-stl-containers"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Conversion from STL containers</h3>
<p>Any sequence container (<code>std::array</code>, <code>std::vector</code>, <code>std::deque</code>, <code>std::forward_list</code>, <code>std::list</code>) whose values can be used to construct JSON values (e.g., integers, floating point numbers, Booleans, string types, or again STL containers described in this section) can be used to create a JSON array. The same holds for similar associative containers (<code>std::set</code>, <code>std::multiset</code>, <code>std::unordered_set</code>, <code>std::unordered_multiset</code>), but in these cases the order of the elements of the array depends on how the elements are ordered in the respective STL container.</p>
<div class="highlight highlight-source-c++"><pre>std::vector&lt;<span class="pl-k">int</span>&gt; c_vector {<span class="pl-c1">1</span>, <span class="pl-c1">2</span>, <span class="pl-c1">3</span>, <span class="pl-c1">4</span>};
json <span class="pl-en">j_vec</span>(c_vector);
<span class="pl-c"><span class="pl-c">//</span> [1, 2, 3, 4]</span>

std::deque&lt;<span class="pl-k">double</span>&gt; c_deque {<span class="pl-c1">1.2</span>, <span class="pl-c1">2.3</span>, <span class="pl-c1">3.4</span>, <span class="pl-c1">5.6</span>};
json <span class="pl-en">j_deque</span>(c_deque);
<span class="pl-c"><span class="pl-c">//</span> [1.2, 2.3, 3.4, 5.6]</span>

std::list&lt;<span class="pl-k">bool</span>&gt; c_list {<span class="pl-c1">true</span>, <span class="pl-c1">true</span>, <span class="pl-c1">false</span>, <span class="pl-c1">true</span>};
json <span class="pl-en">j_list</span>(c_list);
<span class="pl-c"><span class="pl-c">//</span> [true, true, false, true]</span>

std::forward_list&lt;<span class="pl-c1">int64_t</span>&gt; c_flist {<span class="pl-c1">12345678909876</span>, <span class="pl-c1">23456789098765</span>, <span class="pl-c1">34567890987654</span>, <span class="pl-c1">45678909876543</span>};
json <span class="pl-en">j_flist</span>(c_flist);
<span class="pl-c"><span class="pl-c">//</span> [12345678909876, 23456789098765, 34567890987654, 45678909876543]</span>

std::array&lt;<span class="pl-k">unsigned</span> <span class="pl-k">long</span>, <span class="pl-c1">4</span>&gt; c_array {{<span class="pl-c1">1</span>, <span class="pl-c1">2</span>, <span class="pl-c1">3</span>, <span class="pl-c1">4</span>}};
json <span class="pl-en">j_array</span>(c_array);
<span class="pl-c"><span class="pl-c">//</span> [1, 2, 3, 4]</span>

std::set&lt;std::string&gt; c_set {<span class="pl-s"><span class="pl-pds">"</span>one<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>two<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>three<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>four<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>one<span class="pl-pds">"</span></span>};
json <span class="pl-en">j_set</span>(c_set); <span class="pl-c"><span class="pl-c">//</span> only one entry for "one" is used</span>
<span class="pl-c"><span class="pl-c">//</span> ["four", "one", "three", "two"]</span>

std::unordered_set&lt;std::string&gt; c_uset {<span class="pl-s"><span class="pl-pds">"</span>one<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>two<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>three<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>four<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>one<span class="pl-pds">"</span></span>};
json <span class="pl-en">j_uset</span>(c_uset); <span class="pl-c"><span class="pl-c">//</span> only one entry for "one" is used</span>
<span class="pl-c"><span class="pl-c">//</span> maybe ["two", "three", "four", "one"]</span>

std::multiset&lt;std::string&gt; c_mset {<span class="pl-s"><span class="pl-pds">"</span>one<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>two<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>one<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>four<span class="pl-pds">"</span></span>};
json <span class="pl-en">j_mset</span>(c_mset); <span class="pl-c"><span class="pl-c">//</span> both entries for "one" are used</span>
<span class="pl-c"><span class="pl-c">//</span> maybe ["one", "two", "one", "four"]</span>

std::unordered_multiset&lt;std::string&gt; c_umset {<span class="pl-s"><span class="pl-pds">"</span>one<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>two<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>one<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>four<span class="pl-pds">"</span></span>};
json <span class="pl-en">j_umset</span>(c_umset); <span class="pl-c"><span class="pl-c">//</span> both entries for "one" are used</span>
<span class="pl-c"><span class="pl-c">//</span> maybe ["one", "two", "one", "four"]</span></pre></div>
<p>Likewise, any associative key-value containers (<code>std::map</code>, <code>std::multimap</code>, <code>std::unordered_map</code>, <code>std::unordered_multimap</code>) whose keys can construct an <code>std::string</code> and whose values can be used to construct JSON values (see examples above) can be used to create a JSON object. Note that in case of multimaps only one key is used in the JSON object and the value depends on the internal order of the STL container.</p>
<div class="highlight highlight-source-c++"><pre>std::map&lt;std::string, <span class="pl-k">int</span>&gt; c_map { {<span class="pl-s"><span class="pl-pds">"</span>one<span class="pl-pds">"</span></span>, <span class="pl-c1">1</span>}, {<span class="pl-s"><span class="pl-pds">"</span>two<span class="pl-pds">"</span></span>, <span class="pl-c1">2</span>}, {<span class="pl-s"><span class="pl-pds">"</span>three<span class="pl-pds">"</span></span>, <span class="pl-c1">3</span>} };
json <span class="pl-en">j_map</span>(c_map);
<span class="pl-c"><span class="pl-c">//</span> {"one": 1, "three": 3, "two": 2 }</span>

std::unordered_map&lt;<span class="pl-k">const</span> <span class="pl-k">char</span>*, <span class="pl-k">double</span>&gt; c_umap { {<span class="pl-s"><span class="pl-pds">"</span>one<span class="pl-pds">"</span></span>, <span class="pl-c1">1.2</span>}, {<span class="pl-s"><span class="pl-pds">"</span>two<span class="pl-pds">"</span></span>, <span class="pl-c1">2.3</span>}, {<span class="pl-s"><span class="pl-pds">"</span>three<span class="pl-pds">"</span></span>, <span class="pl-c1">3.4</span>} };
json <span class="pl-en">j_umap</span>(c_umap);
<span class="pl-c"><span class="pl-c">//</span> {"one": 1.2, "two": 2.3, "three": 3.4}</span>

std::multimap&lt;std::string, <span class="pl-k">bool</span>&gt; c_mmap { {<span class="pl-s"><span class="pl-pds">"</span>one<span class="pl-pds">"</span></span>, <span class="pl-c1">true</span>}, {<span class="pl-s"><span class="pl-pds">"</span>two<span class="pl-pds">"</span></span>, <span class="pl-c1">true</span>}, {<span class="pl-s"><span class="pl-pds">"</span>three<span class="pl-pds">"</span></span>, <span class="pl-c1">false</span>}, {<span class="pl-s"><span class="pl-pds">"</span>three<span class="pl-pds">"</span></span>, <span class="pl-c1">true</span>} };
json <span class="pl-en">j_mmap</span>(c_mmap); <span class="pl-c"><span class="pl-c">//</span> only one entry for key "three" is used</span>
<span class="pl-c"><span class="pl-c">//</span> maybe {"one": true, "two": true, "three": true}</span>

std::unordered_multimap&lt;std::string, <span class="pl-k">bool</span>&gt; c_ummap { {<span class="pl-s"><span class="pl-pds">"</span>one<span class="pl-pds">"</span></span>, <span class="pl-c1">true</span>}, {<span class="pl-s"><span class="pl-pds">"</span>two<span class="pl-pds">"</span></span>, <span class="pl-c1">true</span>}, {<span class="pl-s"><span class="pl-pds">"</span>three<span class="pl-pds">"</span></span>, <span class="pl-c1">false</span>}, {<span class="pl-s"><span class="pl-pds">"</span>three<span class="pl-pds">"</span></span>, <span class="pl-c1">true</span>} };
json <span class="pl-en">j_ummap</span>(c_ummap); <span class="pl-c"><span class="pl-c">//</span> only one entry for key "three" is used</span>
<span class="pl-c"><span class="pl-c">//</span> maybe {"one": true, "two": true, "three": true}</span></pre></div>
<h3><a id="user-content-json-pointer-and-json-patch" class="anchor" aria-hidden="true" href="#json-pointer-and-json-patch"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>JSON Pointer and JSON Patch</h3>
<p>The library supports <strong>JSON Pointer</strong> (<a href="https://tools.ietf.org/html/rfc6901" rel="nofollow">RFC 6901</a>) as alternative means to address structured values. On top of this, <strong>JSON Patch</strong> (<a href="https://tools.ietf.org/html/rfc6902" rel="nofollow">RFC 6902</a>) allows to describe differences between two JSON values - effectively allowing patch and diff operations known from Unix.</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> a JSON value</span>
json j_original = <span class="pl-s"><span class="pl-pds">R"(</span>{</span>
<span class="pl-s">  "baz": ["one", "two", "three"],</span>
<span class="pl-s">  "foo": "bar"</span>
<span class="pl-s">}<span class="pl-pds">)"</span></span>_json;

<span class="pl-c"><span class="pl-c">//</span> access members with a JSON pointer (RFC 6901)</span>
j_original[<span class="pl-s"><span class="pl-pds">"</span>/baz/1<span class="pl-pds">"</span></span>_json_pointer];
<span class="pl-c"><span class="pl-c">//</span> "two"</span>

<span class="pl-c"><span class="pl-c">//</span> a JSON patch (RFC 6902)</span>
json j_patch = <span class="pl-s"><span class="pl-pds">R"(</span>[</span>
<span class="pl-s">  { "op": "replace", "path": "/baz", "value": "boo" },</span>
<span class="pl-s">  { "op": "add", "path": "/hello", "value": ["world"] },</span>
<span class="pl-s">  { "op": "remove", "path": "/foo"}</span>
<span class="pl-s">]<span class="pl-pds">)"</span></span>_json;

<span class="pl-c"><span class="pl-c">//</span> apply the patch</span>
json j_result = j_original.patch(j_patch);
<span class="pl-c"><span class="pl-c">//</span> {</span>
<span class="pl-c"><span class="pl-c">//</span>    "baz": "boo",</span>
<span class="pl-c"><span class="pl-c">//</span>    "hello": ["world"]</span>
<span class="pl-c"><span class="pl-c">//</span> }</span>

<span class="pl-c"><span class="pl-c">//</span> calculate a JSON patch from two JSON values</span>
<span class="pl-en">json::diff</span>(j_result, j_original);
<span class="pl-c"><span class="pl-c">//</span> [</span>
<span class="pl-c"><span class="pl-c">//</span>   { "op":" replace", "path": "/baz", "value": ["one", "two", "three"] },</span>
<span class="pl-c"><span class="pl-c">//</span>   { "op": "remove","path": "/hello" },</span>
<span class="pl-c"><span class="pl-c">//</span>   { "op": "add", "path": "/foo", "value": "bar" }</span>
<span class="pl-c"><span class="pl-c">//</span> ]</span></pre></div>
<h3><a id="user-content-json-merge-patch" class="anchor" aria-hidden="true" href="#json-merge-patch"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>JSON Merge Patch</h3>
<p>The library supports <strong>JSON Merge Patch</strong> (<a href="https://tools.ietf.org/html/rfc7386" rel="nofollow">RFC 7386</a>) as a patch format. Instead of using JSON Pointer (see above) to specify values to be manipulated, it describes the changes using a syntax that closely mimics the document being modified.</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> a JSON value</span>
json j_document = <span class="pl-s"><span class="pl-pds">R"(</span>{</span>
<span class="pl-s">  "a": "b",</span>
<span class="pl-s">  "c": {</span>
<span class="pl-s">    "d": "e",</span>
<span class="pl-s">    "f": "g"</span>
<span class="pl-s">  }</span>
<span class="pl-s">}<span class="pl-pds">)"</span></span>_json;

<span class="pl-c"><span class="pl-c">//</span> a patch</span>
json j_patch = <span class="pl-s"><span class="pl-pds">R"(</span>{</span>
<span class="pl-s">  "a":"z",</span>
<span class="pl-s">  "c": {</span>
<span class="pl-s">    "f": null</span>
<span class="pl-s">  }</span>
<span class="pl-s">}<span class="pl-pds">)"</span></span>_json;

<span class="pl-c"><span class="pl-c">//</span> apply the patch</span>
j_original.merge_patch(j_patch);
<span class="pl-c"><span class="pl-c">//</span> {</span>
<span class="pl-c"><span class="pl-c">//</span>  "a": "z",</span>
<span class="pl-c"><span class="pl-c">//</span>  "c": {</span>
<span class="pl-c"><span class="pl-c">//</span>    "d": "e"</span>
<span class="pl-c"><span class="pl-c">//</span>  }</span>
<span class="pl-c"><span class="pl-c">//</span> }</span></pre></div>
<h3><a id="user-content-implicit-conversions" class="anchor" aria-hidden="true" href="#implicit-conversions"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Implicit conversions</h3>
<p>The type of the JSON object is determined automatically by the expression to store. Likewise, the stored value is implicitly converted.</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> strings</span>
std::string s1 = <span class="pl-s"><span class="pl-pds">"</span>Hello, world!<span class="pl-pds">"</span></span>;
json js = s1;
std::string s2 = js;

<span class="pl-c"><span class="pl-c">//</span> Booleans</span>
<span class="pl-k">bool</span> b1 = <span class="pl-c1">true</span>;
json jb = b1;
<span class="pl-k">bool</span> b2 = jb;

<span class="pl-c"><span class="pl-c">//</span> numbers</span>
<span class="pl-k">int</span> i = <span class="pl-c1">42</span>;
json jn = i;
<span class="pl-k">double</span> f = jn;

<span class="pl-c"><span class="pl-c">//</span> etc.</span></pre></div>
<p>You can also explicitly ask for the value:</p>
<div class="highlight highlight-source-c++"><pre>std::string vs = js.get&lt;std::string&gt;();
<span class="pl-k">bool</span> vb = jb.get&lt;<span class="pl-k">bool</span>&gt;();
<span class="pl-k">int</span> vi = jn.get&lt;<span class="pl-k">int</span>&gt;();

<span class="pl-c"><span class="pl-c">//</span> etc.</span></pre></div>
<p>Note that <code>char</code> types are not automatically converted to JSON strings, but to integer numbers. A conversion to a string must be specified explicitly:</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-k">char</span> ch = <span class="pl-s"><span class="pl-pds">'</span>A<span class="pl-pds">'</span></span>;                       <span class="pl-c"><span class="pl-c">//</span> ASCII value 65</span>
json j_default = ch;                 <span class="pl-c"><span class="pl-c">//</span> stores integer number 65</span>
json j_string = std::string(<span class="pl-c1">1</span>, ch);  <span class="pl-c"><span class="pl-c">//</span> stores string "A"</span></pre></div>
<h3><a id="user-content-arbitrary-types-conversions" class="anchor" aria-hidden="true" href="#arbitrary-types-conversions"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Arbitrary types conversions</h3>
<p>Every type can be serialized in JSON, not just STL containers and scalar types. Usually, you would do something along those lines:</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-k">namespace</span> <span class="pl-en">ns</span> {
    <span class="pl-c"><span class="pl-c">//</span> a simple struct to model a person</span>
    <span class="pl-k">struct</span> <span class="pl-en">person</span> {
        std::string name;
        std::string address;
        <span class="pl-k">int</span> age;
    };
}

ns::person p = {<span class="pl-s"><span class="pl-pds">"</span>Ned Flanders<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>744 Evergreen Terrace<span class="pl-pds">"</span></span>, <span class="pl-c1">60</span>};

<span class="pl-c"><span class="pl-c">//</span> convert to JSON: copy each value into the JSON object</span>
json j;
j[<span class="pl-s"><span class="pl-pds">"</span>name<span class="pl-pds">"</span></span>] = p.name;
j[<span class="pl-s"><span class="pl-pds">"</span>address<span class="pl-pds">"</span></span>] = p.address;
j[<span class="pl-s"><span class="pl-pds">"</span>age<span class="pl-pds">"</span></span>] = p.age;

<span class="pl-c"><span class="pl-c">//</span> ...</span>

<span class="pl-c"><span class="pl-c">//</span> convert from JSON: copy each value from the JSON object</span>
ns::person p {
    j[<span class="pl-s"><span class="pl-pds">"</span>name<span class="pl-pds">"</span></span>].<span class="pl-smi">get</span>&lt;std::string&gt;(),
    j[<span class="pl-s"><span class="pl-pds">"</span>address<span class="pl-pds">"</span></span>].<span class="pl-smi">get</span>&lt;std::string&gt;(),
    j[<span class="pl-s"><span class="pl-pds">"</span>age<span class="pl-pds">"</span></span>].<span class="pl-smi">get</span>&lt;<span class="pl-k">int</span>&gt;()
};</pre></div>
<p>It works, but that's quite a lot of boilerplate... Fortunately, there's a better way:</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> create a person</span>
ns::person p {<span class="pl-s"><span class="pl-pds">"</span>Ned Flanders<span class="pl-pds">"</span></span>, <span class="pl-s"><span class="pl-pds">"</span>744 Evergreen Terrace<span class="pl-pds">"</span></span>, <span class="pl-c1">60</span>};

<span class="pl-c"><span class="pl-c">//</span> conversion: person -&gt; json</span>
json j = p;

std::cout &lt;&lt; j &lt;&lt; std::endl;
<span class="pl-c"><span class="pl-c">//</span> {"address":"744 Evergreen Terrace","age":60,"name":"Ned Flanders"}</span>

<span class="pl-c"><span class="pl-c">//</span> conversion: json -&gt; person</span>
ns::person p2 = j;

<span class="pl-c"><span class="pl-c">//</span> that's it</span>
<span class="pl-en">assert</span>(p == p2);</pre></div>
<h4><a id="user-content-basic-usage" class="anchor" aria-hidden="true" href="#basic-usage"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Basic usage</h4>
<p>To make this work with one of your types, you only need to provide two functions:</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-k">using</span> nlohmann::json;

<span class="pl-k">namespace</span> <span class="pl-en">ns</span> {
    <span class="pl-k">void</span> <span class="pl-en">to_json</span>(json&amp; j, <span class="pl-k">const</span> person&amp; p) {
        j = json{{<span class="pl-s"><span class="pl-pds">"</span>name<span class="pl-pds">"</span></span>, p.<span class="pl-smi">name</span>}, {<span class="pl-s"><span class="pl-pds">"</span>address<span class="pl-pds">"</span></span>, p.<span class="pl-smi">address</span>}, {<span class="pl-s"><span class="pl-pds">"</span>age<span class="pl-pds">"</span></span>, p.<span class="pl-smi">age</span>}};
    }

    <span class="pl-k">void</span> <span class="pl-en">from_json</span>(<span class="pl-k">const</span> json&amp; j, person&amp; p) {
        j.<span class="pl-c1">at</span>(<span class="pl-s"><span class="pl-pds">"</span>name<span class="pl-pds">"</span></span>).<span class="pl-c1">get_to</span>(p.<span class="pl-smi">name</span>);
        j.<span class="pl-c1">at</span>(<span class="pl-s"><span class="pl-pds">"</span>address<span class="pl-pds">"</span></span>).<span class="pl-c1">get_to</span>(p.<span class="pl-smi">address</span>);
        j.<span class="pl-c1">at</span>(<span class="pl-s"><span class="pl-pds">"</span>age<span class="pl-pds">"</span></span>).<span class="pl-c1">get_to</span>(p.<span class="pl-smi">age</span>);
    }
} <span class="pl-c"><span class="pl-c">//</span> namespace ns</span></pre></div>
<p>That's all! When calling the <code>json</code> constructor with your type, your custom <code>to_json</code> method will be automatically called.
Likewise, when calling <code>get&lt;your_type&gt;()</code> or <code>get_to(your_type&amp;)</code>, the <code>from_json</code> method will be called.</p>
<p>Some important things:</p>
<ul>
<li>Those methods <strong>MUST</strong> be in your type's namespace (which can be the global namespace), or the library will not be able to locate them (in this example, they are in namespace <code>ns</code>, where <code>person</code> is defined).</li>
<li>Those methods <strong>MUST</strong> be available (e.g., properly headers must be included) everywhere you use the implicit conversions. Look at <a href="https://github.com/nlohmann/json/issues/1108">issue 1108</a> for errors that may occur otherwise.</li>
<li>When using <code>get&lt;your_type&gt;()</code>, <code>your_type</code> <strong>MUST</strong> be <a href="https://en.cppreference.com/w/cpp/named_req/DefaultConstructible" rel="nofollow">DefaultConstructible</a>. (There is a way to bypass this requirement described later.)</li>
<li>In function <code>from_json</code>, use function <a href="https://nlohmann.github.io/json/classnlohmann_1_1basic__json_a93403e803947b86f4da2d1fb3345cf2c.html#a93403e803947b86f4da2d1fb3345cf2c" rel="nofollow"><code>at()</code></a> to access the object values rather than <code>operator[]</code>. In case a key does not exist, <code>at</code> throws an exception that you can handle, whereas <code>operator[]</code> exhibits undefined behavior.</li>
<li>In case your type contains several <code>operator=</code> definitions, code like <code>your_variable = your_json;</code> <a href="https://github.com/nlohmann/json/issues/667">may not compile</a>. You need to write <code>your_variable = your_json.get&lt;decltype(your_variable)&gt;();</code> or <code>your_json.get_to(your_variable);</code> instead.</li>
<li>You do not need to add serializers or deserializers for STL types like <code>std::vector</code>: the library already implements these.</li>
</ul>
<h4><a id="user-content-how-do-i-convert-third-party-types" class="anchor" aria-hidden="true" href="#how-do-i-convert-third-party-types"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>How do I convert third-party types?</h4>
<p>This requires a bit more advanced technique. But first, let's see how this conversion mechanism works:</p>
<p>The library uses <strong>JSON Serializers</strong> to convert types to json.
The default serializer for <code>nlohmann::json</code> is <code>nlohmann::adl_serializer</code> (ADL means <a href="https://en.cppreference.com/w/cpp/language/adl" rel="nofollow">Argument-Dependent Lookup</a>).</p>
<p>It is implemented like this (simplified):</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-k">template </span>&lt;<span class="pl-k">typename</span> T&gt;
<span class="pl-k">struct</span> <span class="pl-en">adl_serializer</span> {
    <span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">to_json</span>(json&amp; j, <span class="pl-k">const</span> T&amp; value) {
        <span class="pl-c"><span class="pl-c">//</span> calls the "to_json" method in T's namespace</span>
    }

    <span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">from_json</span>(<span class="pl-k">const</span> json&amp; j, T&amp; value) {
        <span class="pl-c"><span class="pl-c">//</span> same thing, but with the "from_json" method</span>
    }
};</pre></div>
<p>This serializer works fine when you have control over the type's namespace. However, what about <code>boost::optional</code> or <code>std::filesystem::path</code> (C++17)? Hijacking the <code>boost</code> namespace is pretty bad, and it's illegal to add something other than template specializations to <code>std</code>...</p>
<p>To solve this, you need to add a specialization of <code>adl_serializer</code> to the <code>nlohmann</code> namespace, here's an example:</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> partial specialization (full specialization works too)</span>
<span class="pl-k">namespace</span> <span class="pl-en">nlohmann</span> {
    <span class="pl-k">template </span>&lt;<span class="pl-k">typename</span> T&gt;
    <span class="pl-k">struct</span> <span class="pl-en">adl_serializer</span>&lt;boost::optional&lt;T&gt;&gt; {
        <span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">to_json</span>(json&amp; j, <span class="pl-k">const</span> boost::optional&lt;T&gt;&amp; opt) {
            <span class="pl-k">if</span> (opt == boost::none) {
                j = <span class="pl-c1">nullptr</span>;
            } <span class="pl-k">else</span> {
              j = *opt; <span class="pl-c"><span class="pl-c">//</span> this will call adl_serializer&lt;T&gt;::to_json which will</span>
                        <span class="pl-c"><span class="pl-c">//</span> find the free function to_json in T's namespace!</span>
            }
        }

        <span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">from_json</span>(<span class="pl-k">const</span> json&amp; j, boost::optional&lt;T&gt;&amp; opt) {
            <span class="pl-k">if</span> (j.<span class="pl-c1">is_null</span>()) {
                opt = boost::none;
            } <span class="pl-k">else</span> {
                opt = j.<span class="pl-smi">get</span>&lt;T&gt;(); <span class="pl-c"><span class="pl-c">//</span> same as above, but with</span>
                                  <span class="pl-c"><span class="pl-c">//</span> adl_serializer&lt;T&gt;::from_json</span>
            }
        }
    };
}</pre></div>
<h4><a id="user-content-how-can-i-use-get-for-non-default-constructiblenon-copyable-types" class="anchor" aria-hidden="true" href="#how-can-i-use-get-for-non-default-constructiblenon-copyable-types"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>How can I use <code>get()</code> for non-default constructible/non-copyable types?</h4>
<p>There is a way, if your type is <a href="https://en.cppreference.com/w/cpp/named_req/MoveConstructible" rel="nofollow">MoveConstructible</a>. You will need to specialize the <code>adl_serializer</code> as well, but with a special <code>from_json</code> overload:</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-k">struct</span> <span class="pl-en">move_only_type</span> {
    <span class="pl-en">move_only_type</span>() = <span class="pl-k">delete</span>;
    <span class="pl-en">move_only_type</span>(<span class="pl-k">int</span> ii): i(ii) {}
    <span class="pl-en">move_only_type</span>(<span class="pl-k">const</span> move_only_type&amp;) = <span class="pl-k">delete</span>;
    <span class="pl-en">move_only_type</span>(move_only_type&amp;&amp;) = <span class="pl-k">default</span>;

    <span class="pl-k">int</span> i;
};

<span class="pl-k">namespace</span> <span class="pl-en">nlohmann</span> {
    <span class="pl-k">template </span>&lt;&gt;
    <span class="pl-k">struct</span> <span class="pl-en">adl_serializer</span>&lt;move_only_type&gt; {
        <span class="pl-c"><span class="pl-c">//</span> note: the return type is no longer 'void', and the method only takes</span>
        <span class="pl-c"><span class="pl-c">//</span> one argument</span>
        <span class="pl-k">static</span> move_only_type <span class="pl-en">from_json</span>(<span class="pl-k">const</span> json&amp; j) {
            <span class="pl-k">return</span> {j.<span class="pl-smi">get</span>&lt;<span class="pl-k">int</span>&gt;()};
        }

        <span class="pl-c"><span class="pl-c">//</span> Here's the catch! You must provide a to_json method! Otherwise you</span>
        <span class="pl-c"><span class="pl-c">//</span> will not be able to convert move_only_type to json, since you fully</span>
        <span class="pl-c"><span class="pl-c">//</span> specialized adl_serializer on that type</span>
        <span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">to_json</span>(json&amp; j, move_only_type t) {
            j = t.<span class="pl-smi">i</span>;
        }
    };
}</pre></div>
<h4><a id="user-content-can-i-write-my-own-serializer-advanced-use" class="anchor" aria-hidden="true" href="#can-i-write-my-own-serializer-advanced-use"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Can I write my own serializer? (Advanced use)</h4>
<p>Yes. You might want to take a look at <a href="https://github.com/nlohmann/json/blob/develop/test/src/unit-udt.cpp"><code>unit-udt.cpp</code></a> in the test suite, to see a few examples.</p>
<p>If you write your own serializer, you'll need to do a few things:</p>
<ul>
<li>use a different <code>basic_json</code> alias than <code>nlohmann::json</code> (the last template parameter of <code>basic_json</code> is the <code>JSONSerializer</code>)</li>
<li>use your <code>basic_json</code> alias (or a template parameter) in all your <code>to_json</code>/<code>from_json</code> methods</li>
<li>use <code>nlohmann::to_json</code> and <code>nlohmann::from_json</code> when you need ADL</li>
</ul>
<p>Here is an example, without simplifications, that only accepts types with a size &lt;= 32, and uses ADL.</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> You should use void as a second template argument</span>
<span class="pl-c"><span class="pl-c">//</span> if you don't need compile-time checks on T</span>
<span class="pl-k">template</span>&lt;<span class="pl-k">typename</span> T, <span class="pl-k">typename</span> SFINAE = <span class="pl-k">typename</span> std::enable_if&lt;<span class="pl-k">sizeof</span>(T) &lt;= <span class="pl-c1">32</span>&gt;::type&gt;
<span class="pl-k">struct</span> <span class="pl-en">less_than_32_serializer</span> {
    <span class="pl-k">template </span>&lt;<span class="pl-k">typename</span> BasicJsonType&gt;
    <span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">to_json</span>(BasicJsonType&amp; j, T value) {
        <span class="pl-c"><span class="pl-c">//</span> we want to use ADL, and call the correct to_json overload</span>
        <span class="pl-k">using</span> nlohmann::to_json; <span class="pl-c"><span class="pl-c">//</span> this method is called by adl_serializer,</span>
                                 <span class="pl-c"><span class="pl-c">//</span> this is where the magic happens</span>
        <span class="pl-c1">to_json</span>(j, value);
    }

    <span class="pl-k">template </span>&lt;<span class="pl-k">typename</span> BasicJsonType&gt;
    <span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">from_json</span>(<span class="pl-k">const</span> BasicJsonType&amp; j, T&amp; value) {
        <span class="pl-c"><span class="pl-c">//</span> same thing here</span>
        <span class="pl-k">using</span> nlohmann::from_json;
        <span class="pl-c1">from_json</span>(j, value);
    }
};</pre></div>
<p>Be <strong>very</strong> careful when reimplementing your serializer, you can stack overflow if you don't pay attention:</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-k">template </span>&lt;<span class="pl-k">typename</span> T, <span class="pl-k">void</span>&gt;
<span class="pl-k">struct</span> <span class="pl-en">bad_serializer</span>
{
    <span class="pl-k">template </span>&lt;<span class="pl-k">typename</span> BasicJsonType&gt;
    <span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">to_json</span>(BasicJsonType&amp; j, <span class="pl-k">const</span> T&amp; value) {
      <span class="pl-c"><span class="pl-c">//</span> this calls BasicJsonType::json_serializer&lt;T&gt;::to_json(j, value);</span>
      <span class="pl-c"><span class="pl-c">//</span> if BasicJsonType::json_serializer == bad_serializer ... oops!</span>
      j = value;
    }

    <span class="pl-k">template </span>&lt;<span class="pl-k">typename</span> BasicJsonType&gt;
    <span class="pl-k">static</span> <span class="pl-k">void</span> <span class="pl-en">to_json</span>(<span class="pl-k">const</span> BasicJsonType&amp; j, T&amp; value) {
      <span class="pl-c"><span class="pl-c">//</span> this calls BasicJsonType::json_serializer&lt;T&gt;::from_json(j, value);</span>
      <span class="pl-c"><span class="pl-c">//</span> if BasicJsonType::json_serializer == bad_serializer ... oops!</span>
      value = j.<span class="pl-k">template</span> <span class="pl-smi">get</span>&lt;T&gt;(); <span class="pl-c"><span class="pl-c">//</span> oops!</span>
    }
};</pre></div>
<h3><a id="user-content-specializing-enum-conversion" class="anchor" aria-hidden="true" href="#specializing-enum-conversion"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Specializing enum conversion</h3>
<p>By default, enum values are serialized to JSON as integers. In some cases this could result in undesired behavior. If an enum is modified or re-ordered after data has been serialized to JSON, the later de-serialized JSON data may be undefined or a different enum value than was originally intended.</p>
<p>It is possible to more precisely specify how a given enum is mapped to and from JSON as shown below:</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> example enum type declaration</span>
<span class="pl-k">enum</span> TaskState {
    TS_STOPPED,
    TS_RUNNING,
    TS_COMPLETED,
    TS_INVALID=-<span class="pl-c1">1</span>,
};

<span class="pl-c"><span class="pl-c">//</span> map TaskState values to JSON as strings</span>
<span class="pl-en">NLOHMANN_JSON_SERIALIZE_ENUM</span>( TaskState, {
    {TS_INVALID, <span class="pl-c1">nullptr</span>},
    {TS_STOPPED, <span class="pl-s"><span class="pl-pds">"</span>stopped<span class="pl-pds">"</span></span>},
    {TS_RUNNING, <span class="pl-s"><span class="pl-pds">"</span>running<span class="pl-pds">"</span></span>},
    {TS_COMPLETED, <span class="pl-s"><span class="pl-pds">"</span>completed<span class="pl-pds">"</span></span>},
});</pre></div>
<p>The <code>NLOHMANN_JSON_SERIALIZE_ENUM()</code> macro declares a set of <code>to_json()</code> / <code>from_json()</code> functions for type <code>TaskState</code> while avoiding repetition and boilerplate serilization code.</p>
<p><strong>Usage:</strong></p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> enum to JSON as string</span>
json j = TS_STOPPED;
<span class="pl-en">assert</span>(j == <span class="pl-s"><span class="pl-pds">"</span>stopped<span class="pl-pds">"</span></span>);

<span class="pl-c"><span class="pl-c">//</span> json string to enum</span>
json j3 = <span class="pl-s"><span class="pl-pds">"</span>running<span class="pl-pds">"</span></span>;
<span class="pl-en">assert</span>(j3.get&lt;TaskState&gt;() == TS_RUNNING);

<span class="pl-c"><span class="pl-c">//</span> undefined json value to enum (where the first map entry above is the default)</span>
json jPi = <span class="pl-c1">3.14</span>;
<span class="pl-en">assert</span>(jPi.get&lt;TaskState&gt;() == TS_INVALID );</pre></div>
<p>Just as in <a href="#arbitrary-types-conversions">Arbitrary Type Conversions</a> above,</p>
<ul>
<li><code>NLOHMANN_JSON_SERIALIZE_ENUM()</code> MUST be declared in your enum type's namespace (which can be the global namespace), or the library will not be able to locate it and it will default to integer serialization.</li>
<li>It MUST be available (e.g., proper headers must be included) everywhere you use the conversions.</li>
</ul>
<p>Other Important points:</p>
<ul>
<li>When using <code>get&lt;ENUM_TYPE&gt;()</code>, undefined JSON values will default to the first pair specified in your map. Select this default pair carefully.</li>
<li>If an enum or JSON value is specified more than once in your map, the first matching occurrence from the top of the map will be returned when converting to or from JSON.</li>
</ul>
<h3><a id="user-content-binary-formats-bson-cbor-messagepack-and-ubjson" class="anchor" aria-hidden="true" href="#binary-formats-bson-cbor-messagepack-and-ubjson"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Binary formats (BSON, CBOR, MessagePack, and UBJSON</h3>
<p>Though JSON is a ubiquitous data format, it is not a very compact format suitable for data exchange, for instance over a network. Hence, the library supports <a href="http://bsonspec.org" rel="nofollow">BSON</a> (Binary JSON), <a href="http://cbor.io" rel="nofollow">CBOR</a> (Concise Binary Object Representation), <a href="http://msgpack.org" rel="nofollow">MessagePack</a>, and <a href="http://ubjson.org" rel="nofollow">UBJSON</a> (Universal Binary JSON Specification) to efficiently encode JSON values to byte vectors and to decode such vectors.</p>
<div class="highlight highlight-source-c++"><pre><span class="pl-c"><span class="pl-c">//</span> create a JSON value</span>
json j = <span class="pl-s"><span class="pl-pds">R"(</span>{"compact": true, "schema": 0}<span class="pl-pds">)"</span></span>_json;

<span class="pl-c"><span class="pl-c">//</span> serialize to BSON</span>
std::vector&lt;std::<span class="pl-c1">uint8_t</span>&gt; v_bson = json::to_bson(j);

<span class="pl-c"><span class="pl-c">//</span> 0x1B, 0x00, 0x00, 0x00, 0x08, 0x63, 0x6F, 0x6D, 0x70, 0x61, 0x63, 0x74, 0x00, 0x01, 0x10, 0x73, 0x63, 0x68, 0x65, 0x6D, 0x61, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00</span>

<span class="pl-c"><span class="pl-c">//</span> roundtrip</span>
json j_from_bson = json::from_bson(v_bson);

<span class="pl-c"><span class="pl-c">//</span> serialize to CBOR</span>
std::vector&lt;std::<span class="pl-c1">uint8_t</span>&gt; v_cbor = json::to_cbor(j);

<span class="pl-c"><span class="pl-c">//</span> 0xA2, 0x67, 0x63, 0x6F, 0x6D, 0x70, 0x61, 0x63, 0x74, 0xF5, 0x66, 0x73, 0x63, 0x68, 0x65, 0x6D, 0x61, 0x00</span>

<span class="pl-c"><span class="pl-c">//</span> roundtrip</span>
json j_from_cbor = json::from_cbor(v_cbor);

<span class="pl-c"><span class="pl-c">//</span> serialize to MessagePack</span>
std::vector&lt;std::<span class="pl-c1">uint8_t</span>&gt; v_msgpack = json::to_msgpack(j);

<span class="pl-c"><span class="pl-c">//</span> 0x82, 0xA7, 0x63, 0x6F, 0x6D, 0x70, 0x61, 0x63, 0x74, 0xC3, 0xA6, 0x73, 0x63, 0x68, 0x65, 0x6D, 0x61, 0x00</span>

<span class="pl-c"><span class="pl-c">//</span> roundtrip</span>
json j_from_msgpack = json::from_msgpack(v_msgpack);

<span class="pl-c"><span class="pl-c">//</span> serialize to UBJSON</span>
std::vector&lt;std::<span class="pl-c1">uint8_t</span>&gt; v_ubjson = json::to_ubjson(j);

<span class="pl-c"><span class="pl-c">//</span> 0x7B, 0x69, 0x07, 0x63, 0x6F, 0x6D, 0x70, 0x61, 0x63, 0x74, 0x54, 0x69, 0x06, 0x73, 0x63, 0x68, 0x65, 0x6D, 0x61, 0x69, 0x00, 0x7D</span>

<span class="pl-c"><span class="pl-c">//</span> roundtrip</span>
json j_from_ubjson = json::from_ubjson(v_ubjson);</pre></div>
<h2><a id="user-content-supported-compilers" class="anchor" aria-hidden="true" href="#supported-compilers"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Supported compilers</h2>
<p>Though it's 2018 already, the support for C++11 is still a bit sparse. Currently, the following compilers are known to work:</p>
<ul>
<li>GCC 4.8 - 8.2 (and possibly later)</li>
<li>Clang 3.4 - 6.1 (and possibly later)</li>
<li>Intel C++ Compiler 17.0.2 (and possibly later)</li>
<li>Microsoft Visual C++ 2015 / Build Tools 14.0.25123.0 (and possibly later)</li>
<li>Microsoft Visual C++ 2017 / Build Tools 15.5.180.51428 (and possibly later)</li>
</ul>
<p>I would be happy to learn about other compilers/versions.</p>
<p>Please note:</p>
<ul>
<li>
<p>GCC 4.8 has a bug <a href="https://gcc.gnu.org/bugzilla/show_bug.cgi?id=57824" rel="nofollow">57824</a>): multiline raw strings cannot be the arguments to macros. Don't use multiline raw strings directly in macros with this compiler.</p>
</li>
<li>
<p>Android defaults to using very old compilers and C++ libraries. To fix this, add the following to your <code>Application.mk</code>. This will switch to the LLVM C++ library, the Clang compiler, and enable C++11 and other features disabled by default.</p>
<pre><code>APP_STL := c++_shared
NDK_TOOLCHAIN_VERSION := clang3.6
APP_CPPFLAGS += -frtti -fexceptions
</code></pre>
<p>The code compiles successfully with <a href="https://developer.android.com/ndk/index.html?hl=ml" rel="nofollow">Android NDK</a>, Revision 9 - 11 (and possibly later) and <a href="https://www.crystax.net/en/android/ndk" rel="nofollow">CrystaX's Android NDK</a> version 10.</p>
</li>
<li>
<p>For GCC running on MinGW or Android SDK, the error <code>'to_string' is not a member of 'std'</code> (or similarly, for <code>strtod</code>) may occur. Note this is not an issue with the code,  but rather with the compiler itself. On Android, see above to build with a newer environment.  For MinGW, please refer to <a href="http://tehsausage.com/mingw-to-string" rel="nofollow">this site</a> and <a href="https://github.com/nlohmann/json/issues/136">this discussion</a> for information on how to fix this bug. For Android NDK using <code>APP_STL := gnustl_static</code>, please refer to <a href="https://github.com/nlohmann/json/issues/219">this discussion</a>.</p>
</li>
<li>
<p>Unsupported versions of GCC and Clang are rejected by <code>#error</code> directives. This can be switched off by defining <code>JSON_SKIP_UNSUPPORTED_COMPILER_CHECK</code>. Note that you can expect no support in this case.</p>
</li>
</ul>
<p>The following compilers are currently used in continuous integration at <a href="https://travis-ci.org/nlohmann/json" rel="nofollow">Travis</a> and <a href="https://ci.appveyor.com/project/nlohmann/json" rel="nofollow">AppVeyor</a>:</p>
<table>
<thead>
<tr>
<th>Compiler</th>
<th>Operating System</th>
<th>Version String</th>
</tr>
</thead>
<tbody>
<tr>
<td>GCC 4.8.5</td>
<td>Ubuntu 14.04.5 LTS</td>
<td>g++-4.8 (Ubuntu 4.8.5-2ubuntu1~14.04.2) 4.8.5</td>
</tr>
<tr>
<td>GCC 4.9.4</td>
<td>Ubuntu 14.04.1 LTS</td>
<td>g++-4.9 (Ubuntu 4.9.4-2ubuntu1~14.04.1) 4.9.4</td>
</tr>
<tr>
<td>GCC 5.5.0</td>
<td>Ubuntu 14.04.1 LTS</td>
<td>g++-5 (Ubuntu 5.5.0-12ubuntu1~14.04) 5.5.0 20171010</td>
</tr>
<tr>
<td>GCC 6.4.0</td>
<td>Ubuntu 14.04.1 LTS</td>
<td>g++-6 (Ubuntu 6.4.0-17ubuntu1~14.04) 6.4.0 20180424</td>
</tr>
<tr>
<td>GCC 7.3.0</td>
<td>Ubuntu 14.04.1 LTS</td>
<td>g++-7 (Ubuntu 7.3.0-21ubuntu1~14.04) 7.3.0</td>
</tr>
<tr>
<td>GCC 7.3.0</td>
<td>Windows Server 2012 R2 (x64)</td>
<td>g++ (x86_64-posix-seh-rev0, Built by MinGW-W64 project) 7.3.0</td>
</tr>
<tr>
<td>GCC 8.1.0</td>
<td>Ubuntu 14.04.1 LTS</td>
<td>g++-8 (Ubuntu 8.1.0-5ubuntu1~14.04) 8.1.0</td>
</tr>
<tr>
<td>Clang 3.5.0</td>
<td>Ubuntu 14.04.1 LTS</td>
<td>clang version 3.5.0-4ubuntu2~trusty2 (tags/RELEASE_350/final) (based on LLVM 3.5.0)</td>
</tr>
<tr>
<td>Clang 3.6.2</td>
<td>Ubuntu 14.04.1 LTS</td>
<td>clang version 3.6.2-svn240577-1~exp1 (branches/release_36) (based on LLVM 3.6.2)</td>
</tr>
<tr>
<td>Clang 3.7.1</td>
<td>Ubuntu 14.04.1 LTS</td>
<td>clang version 3.7.1-svn253571-1~exp1 (branches/release_37) (based on LLVM 3.7.1)</td>
</tr>
<tr>
<td>Clang 3.8.0</td>
<td>Ubuntu 14.04.1 LTS</td>
<td>clang version 3.8.0-2ubuntu3~trusty5 (tags/RELEASE_380/final)</td>
</tr>
<tr>
<td>Clang 3.9.1</td>
<td>Ubuntu 14.04.1 LTS</td>
<td>clang version 3.9.1-4ubuntu3~14.04.3 (tags/RELEASE_391/rc2)</td>
</tr>
<tr>
<td>Clang 4.0.1</td>
<td>Ubuntu 14.04.1 LTS</td>
<td>clang version 4.0.1-svn305264-1~exp1 (branches/release_40)</td>
</tr>
<tr>
<td>Clang 5.0.2</td>
<td>Ubuntu 14.04.1 LTS</td>
<td>clang version 5.0.2-svn328729-1<del>exp1</del>20180509123505.100 (branches/release_50)</td>
</tr>
<tr>
<td>Clang 6.0.1</td>
<td>Ubuntu 14.04.1 LTS</td>
<td>clang version 6.0.1-svn334776-1<del>exp1</del>20180726133705.85 (branches/release_60)</td>
</tr>
<tr>
<td>Clang Xcode 6.4</td>
<td>OSX 10.10.5</td>
<td>Apple LLVM version 6.1.0 (clang-602.0.53) (based on LLVM 3.6.0svn)</td>
</tr>
<tr>
<td>Clang Xcode 7.3</td>
<td>OSX 10.11.6</td>
<td>Apple LLVM version 7.3.0 (clang-703.0.31)</td>
</tr>
<tr>
<td>Clang Xcode 8.0</td>
<td>OSX 10.11.6</td>
<td>Apple LLVM version 8.0.0 (clang-800.0.38)</td>
</tr>
<tr>
<td>Clang Xcode 8.1</td>
<td>OSX 10.12.6</td>
<td>Apple LLVM version 8.0.0 (clang-800.0.42.1)</td>
</tr>
<tr>
<td>Clang Xcode 8.2</td>
<td>OSX 10.12.6</td>
<td>Apple LLVM version 8.0.0 (clang-800.0.42.1)</td>
</tr>
<tr>
<td>Clang Xcode 8.3</td>
<td>OSX 10.11.6</td>
<td>Apple LLVM version 8.1.0 (clang-802.0.38)</td>
</tr>
<tr>
<td>Clang Xcode 9.0</td>
<td>OSX 10.12.6</td>
<td>Apple LLVM version 9.0.0 (clang-900.0.37)</td>
</tr>
<tr>
<td>Clang Xcode 9.1</td>
<td>OSX 10.12.6</td>
<td>Apple LLVM version 9.0.0 (clang-900.0.38)</td>
</tr>
<tr>
<td>Clang Xcode 9.2</td>
<td>OSX 10.13.3</td>
<td>Apple LLVM version 9.1.0 (clang-902.0.39.1)</td>
</tr>
<tr>
<td>Clang Xcode 9.3</td>
<td>OSX 10.13.3</td>
<td>Apple LLVM version 9.1.0 (clang-902.0.39.2)</td>
</tr>
<tr>
<td>Clang Xcode 10.0</td>
<td>OSX 10.13.3</td>
<td>Apple LLVM version 10.0.0 (clang-1000.11.45.2)</td>
</tr>
<tr>
<td>Visual Studio 14 2015</td>
<td>Windows Server 2012 R2 (x64)</td>
<td>Microsoft (R) Build Engine version 14.0.25420.1, MSVC 19.0.24215.1</td>
</tr>
<tr>
<td>Visual Studio 2017</td>
<td>Windows Server 2016</td>
<td>Microsoft (R) Build Engine version 15.7.180.61344, MSVC 19.14.26433.0</td>
</tr>
</tbody>
</table>
<h2><a id="user-content-license" class="anchor" aria-hidden="true" href="#license"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>License</h2>
<p><a target="_blank" rel="noopener noreferrer" href="https://camo.githubusercontent.com/d9f2a52ccb094aecca865c7614750675ddf80fdb/687474703a2f2f6f70656e736f757263652e6f72672f74726164656d61726b732f6f70656e736f757263652f4f53492d417070726f7665642d4c6963656e73652d313030783133372e706e67"><img align="right" src="https://camo.githubusercontent.com/d9f2a52ccb094aecca865c7614750675ddf80fdb/687474703a2f2f6f70656e736f757263652e6f72672f74726164656d61726b732f6f70656e736f757263652f4f53492d417070726f7665642d4c6963656e73652d313030783133372e706e67" data-canonical-src="http://opensource.org/trademarks/opensource/OSI-Approved-License-100x137.png" style="max-width:100%;"></a></p>
<p>The class is licensed under the <a href="http://opensource.org/licenses/MIT" rel="nofollow">MIT License</a>:</p>
<p>Copyright © 2013-2018 <a href="http://nlohmann.me" rel="nofollow">Niels Lohmann</a></p>
<p>Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:</p>
<p>The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.</p>
<p>THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.</p>
<hr>
<p>The class contains the UTF-8 Decoder from Bjoern Hoehrmann which is licensed under the <a href="http://opensource.org/licenses/MIT" rel="nofollow">MIT License</a> (see above). Copyright © 2008-2009 <a href="http://bjoern.hoehrmann.de/" rel="nofollow">Björn Hoehrmann</a> <a href="mailto:bjoern@hoehrmann.de">bjoern@hoehrmann.de</a></p>
<p>The class contains a slightly modified version of the Grisu2 algorithm from Florian Loitsch which is licensed under the <a href="http://opensource.org/licenses/MIT" rel="nofollow">MIT License</a> (see above). Copyright © 2009 <a href="http://florian.loitsch.com/" rel="nofollow">Florian Loitsch</a></p>
<h2><a id="user-content-contact" class="anchor" aria-hidden="true" href="#contact"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Contact</h2>
<p>If you have questions regarding the library, I would like to invite you to <a href="https://github.com/nlohmann/json/issues/new">open an issue at GitHub</a>. Please describe your request, problem, or question as detailed as possible, and also mention the version of the library you are using as well as the version of your compiler and operating system. Opening an issue at GitHub allows other users and contributors to this library to collaborate. For instance, I have little experience with MSVC, and most issues in this regard have been solved by a growing community. If you have a look at the <a href="https://github.com/nlohmann/json/issues?q=is%3Aissue+is%3Aclosed">closed issues</a>, you will see that we react quite timely in most cases.</p>
<p>Only if your request would contain confidential information, please <a href="mailto:mail@nlohmann.me">send me an email</a>. For encrypted messages, please use <a href="https://keybase.io/nlohmann/pgp_keys.asc" rel="nofollow">this key</a>.</p>
<h2><a id="user-content-security" class="anchor" aria-hidden="true" href="#security"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Security</h2>
<p><a href="https://github.com/nlohmann/json/commits">Commits by Niels Lohmann</a> and <a href="https://github.com/nlohmann/json/releases">releases</a> are signed with this <a href="https://keybase.io/nlohmann/pgp_keys.asc?fingerprint=797167ae41c0a6d9232e48457f3cea63ae251b69" rel="nofollow">PGP Key</a>.</p>
<h2><a id="user-content-thanks" class="anchor" aria-hidden="true" href="#thanks"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Thanks</h2>
<p>I deeply appreciate the help of the following people.</p>
<p><a target="_blank" rel="noopener noreferrer" href="https://raw.githubusercontent.com/nlohmann/json/develop/doc/avatars.png"><img src="https://raw.githubusercontent.com/nlohmann/json/develop/doc/avatars.png" alt="Contributors" style="max-width:100%;"></a></p>
<ul>
<li><a href="https://github.com/Teemperor">Teemperor</a> implemented CMake support and lcov integration, realized escape and Unicode handling in the string parser, and fixed the JSON serialization.</li>
<li><a href="https://github.com/elliotgoodrich">elliotgoodrich</a> fixed an issue with double deletion in the iterator classes.</li>
<li><a href="https://github.com/kirkshoop">kirkshoop</a> made the iterators of the class composable to other libraries.</li>
<li><a href="https://github.com/wanwc">wancw</a> fixed a bug that hindered the class to compile with Clang.</li>
<li>Tomas Åblad found a bug in the iterator implementation.</li>
<li><a href="https://github.com/jrandall">Joshua C. Randall</a> fixed a bug in the floating-point serialization.</li>
<li><a href="https://github.com/aburgh">Aaron Burghardt</a> implemented code to parse streams incrementally. Furthermore, he greatly improved the parser class by allowing the definition of a filter function to discard undesired elements while parsing.</li>
<li><a href="https://github.com/dkopecek">Daniel Kopeček</a> fixed a bug in the compilation with GCC 5.0.</li>
<li><a href="https://github.com/Florianjw">Florian Weber</a> fixed a bug in and improved the performance of the comparison operators.</li>
<li><a href="https://github.com/EricMCornelius">Eric Cornelius</a> pointed out a bug in the handling with NaN and infinity values. He also improved the performance of the string escaping.</li>
<li><a href="https://github.com/likebeta">易思龙</a> implemented a conversion from anonymous enums.</li>
<li><a href="https://github.com/kepkin">kepkin</a> patiently pushed forward the support for Microsoft Visual studio.</li>
<li><a href="https://github.com/gregmarr">gregmarr</a> simplified the implementation of reverse iterators and helped with numerous hints and improvements. In particular, he pushed forward the implementation of user-defined types.</li>
<li><a href="https://github.com/caiovlp">Caio Luppi</a> fixed a bug in the Unicode handling.</li>
<li><a href="https://github.com/dariomt">dariomt</a> fixed some typos in the examples.</li>
<li><a href="https://github.com/d-frey">Daniel Frey</a> cleaned up some pointers and implemented exception-safe memory allocation.</li>
<li><a href="https://github.com/ColinH">Colin Hirsch</a> took care of a small namespace issue.</li>
<li><a href="https://github.com/whoshuu">Huu Nguyen</a> correct a variable name in the documentation.</li>
<li><a href="https://github.com/silverweed">Silverweed</a> overloaded <code>parse()</code> to accept an rvalue reference.</li>
<li><a href="https://github.com/dariomt">dariomt</a> fixed a subtlety in MSVC type support and implemented the <code>get_ref()</code> function to get a reference to stored values.</li>
<li><a href="https://github.com/ZahlGraf">ZahlGraf</a> added a workaround that allows compilation using Android NDK.</li>
<li><a href="https://github.com/whackashoe">whackashoe</a> replaced a function that was marked as unsafe by Visual Studio.</li>
<li><a href="https://github.com/406345">406345</a> fixed two small warnings.</li>
<li><a href="https://github.com/glenfe">Glen Fernandes</a> noted a potential portability problem in the <code>has_mapped_type</code> function.</li>
<li><a href="https://github.com/nibroc">Corbin Hughes</a> fixed some typos in the contribution guidelines.</li>
<li><a href="https://github.com/twelsby">twelsby</a> fixed the array subscript operator, an issue that failed the MSVC build, and floating-point parsing/dumping. He further added support for unsigned integer numbers and implemented better roundtrip support for parsed numbers.</li>
<li><a href="https://github.com/vog">Volker Diels-Grabsch</a> fixed a link in the README file.</li>
<li><a href="https://github.com/msm-">msm-</a> added support for American Fuzzy Lop.</li>
<li><a href="https://github.com/Annihil">Annihil</a> fixed an example in the README file.</li>
<li><a href="https://github.com/Themercee">Themercee</a> noted a wrong URL in the README file.</li>
<li><a href="https://github.com/lv-zheng">Lv Zheng</a> fixed a namespace issue with <code>int64_t</code> and <code>uint64_t</code>.</li>
<li><a href="https://github.com/abc100m">abc100m</a> analyzed the issues with GCC 4.8 and proposed a <a href="https://github.com/nlohmann/json/pull/212">partial solution</a>.</li>
<li><a href="https://github.com/zewt">zewt</a> added useful notes to the README file about Android.</li>
<li><a href="https://github.com/robertmrk">Róbert Márki</a> added a fix to use move iterators and improved the integration via CMake.</li>
<li><a href="https://github.com/ChrisKitching">Chris Kitching</a> cleaned up the CMake files.</li>
<li><a href="https://github.com/06needhamt">Tom Needham</a> fixed a subtle bug with MSVC 2015 which was also proposed by <a href="https://github.com/Epidal">Michael K.</a>.</li>
<li><a href="https://github.com/thelostt">Mário Feroldi</a> fixed a small typo.</li>
<li><a href="https://github.com/duncanwerner">duncanwerner</a> found a really embarrassing performance regression in the 2.0.0 release.</li>
<li><a href="https://github.com/dtoma">Damien</a> fixed one of the last conversion warnings.</li>
<li><a href="https://github.com/t-b">Thomas Braun</a> fixed a warning in a test case.</li>
<li><a href="https://github.com/theodelrieu">Théo DELRIEU</a> patiently and constructively oversaw the long way toward <a href="https://github.com/nlohmann/json/issues/290">iterator-range parsing</a>. He also implemented the magic behind the serialization/deserialization of user-defined types and split the single header file into smaller chunks.</li>
<li><a href="https://github.com/5tefan">Stefan</a> fixed a minor issue in the documentation.</li>
<li><a href="https://github.com/vasild">Vasil Dimov</a> fixed the documentation regarding conversions from <code>std::multiset</code>.</li>
<li><a href="https://github.com/ChristophJud">ChristophJud</a> overworked the CMake files to ease project inclusion.</li>
<li><a href="https://github.com/vpetrigo">Vladimir Petrigo</a> made a SFINAE hack more readable and added Visual Studio 17 to the build matrix.</li>
<li><a href="https://github.com/seeekr">Denis Andrejew</a> fixed a grammar issue in the README file.</li>
<li><a href="https://github.com/palacaze">Pierre-Antoine Lacaze</a> found a subtle bug in the <code>dump()</code> function.</li>
<li><a href="https://github.com/TurpentineDistillery">TurpentineDistillery</a> pointed to <a href="https://en.cppreference.com/w/cpp/locale/locale/classic" rel="nofollow"><code>std::locale::classic()</code></a> to avoid too much locale joggling, found some nice performance improvements in the parser, improved the benchmarking code, and realized locale-independent number parsing and printing.</li>
<li><a href="https://github.com/cgzones">cgzones</a> had an idea how to fix the Coverity scan.</li>
<li><a href="https://github.com/jaredgrubb">Jared Grubb</a> silenced a nasty documentation warning.</li>
<li><a href="https://github.com/qwename">Yixin Zhang</a> fixed an integer overflow check.</li>
<li><a href="https://github.com/Bosswestfalen">Bosswestfalen</a> merged two iterator classes into a smaller one.</li>
<li><a href="https://github.com/Daniel599">Daniel599</a> helped to get Travis execute the tests with Clang's sanitizers.</li>
<li><a href="https://github.com/vjon">Jonathan Lee</a> fixed an example in the README file.</li>
<li><a href="https://github.com/gnzlbg">gnzlbg</a> supported the implementation of user-defined types.</li>
<li><a href="https://github.com/qis">Alexej Harm</a> helped to get the user-defined types working with Visual Studio.</li>
<li><a href="https://github.com/jaredgrubb">Jared Grubb</a> supported the implementation of user-defined types.</li>
<li><a href="https://github.com/EnricoBilla">EnricoBilla</a> noted a typo in an example.</li>
<li><a href="https://github.com/horenmar">Martin Hořeňovský</a> found a way for a 2x speedup for the compilation time of the test suite.</li>
<li><a href="https://github.com/ukhegg">ukhegg</a> found proposed an improvement for the examples section.</li>
<li><a href="https://github.com/rswanson-ihi">rswanson-ihi</a> noted a typo in the README.</li>
<li><a href="https://github.com/stanmihai4">Mihai Stan</a> fixed a bug in the comparison with <code>nullptr</code>s.</li>
<li><a href="https://github.com/tusharpm">Tushar Maheshwari</a> added <a href="https://github.com/sakra/cotire">cotire</a> support to speed up the compilation.</li>
<li><a href="https://github.com/TedLyngmo">TedLyngmo</a> noted a typo in the README, removed unnecessary bit arithmetic, and fixed some <code>-Weffc++</code> warnings.</li>
<li><a href="https://github.com/krzysztofwos">Krzysztof Woś</a> made exceptions more visible.</li>
<li><a href="https://github.com/ftillier">ftillier</a> fixed a compiler warning.</li>
<li><a href="https://github.com/tinloaf">tinloaf</a> made sure all pushed warnings are properly popped.</li>
<li><a href="https://github.com/Fytch">Fytch</a> found a bug in the documentation.</li>
<li><a href="https://github.com/Type1J">Jay Sistar</a> implemented a Meson build description.</li>
<li><a href="https://github.com/HenryRLee">Henry Lee</a> fixed a warning in ICC and improved the iterator implementation.</li>
<li><a href="https://github.com/vthiery">Vincent Thiery</a> maintains a package for the Conan package manager.</li>
<li><a href="https://github.com/koemeet">Steffen</a> fixed a potential issue with MSVC and <code>std::min</code>.</li>
<li><a href="https://github.com/Chocobo1">Mike Tzou</a> fixed some typos.</li>
<li><a href="https://github.com/amrcode">amrcode</a> noted a misleading documentation about comparison of floats.</li>
<li><a href="https://github.com/olegendo">Oleg Endo</a> reduced the memory consumption by replacing <code>&lt;iostream&gt;</code> with <code>&lt;iosfwd&gt;</code>.</li>
<li><a href="https://github.com/dan-42">dan-42</a> cleaned up the CMake files to simplify including/reusing of the library.</li>
<li><a href="https://github.com/himikof">Nikita Ofitserov</a> allowed for moving values from initializer lists.</li>
<li><a href="https://github.com/wincent">Greg Hurrell</a> fixed a typo.</li>
<li><a href="https://github.com/DmitryKuk">Dmitry Kukovinets</a> fixed a typo.</li>
<li><a href="https://github.com/kbthomp1">kbthomp1</a> fixed an issue related to the Intel OSX compiler.</li>
<li><a href="https://github.com/daixtrose">Markus Werle</a> fixed a typo.</li>
<li><a href="https://github.com/WebProdPP">WebProdPP</a> fixed a subtle error in a precondition check.</li>
<li><a href="https://github.com/leha-bot">Alex</a> noted an error in a code sample.</li>
<li><a href="https://github.com/tdegeus">Tom de Geus</a> reported some warnings with ICC and helped fixing them.</li>
<li><a href="https://github.com/pjkundert">Perry Kundert</a> simplified reading from input streams.</li>
<li><a href="https://github.com/sonulohani">Sonu Lohani</a> fixed a small compilation error.</li>
<li><a href="https://github.com/jseward">Jamie Seward</a> fixed all MSVC warnings.</li>
<li><a href="https://github.com/eld00d">Nate Vargas</a> added a Doxygen tag file.</li>
<li><a href="https://github.com/pvleuven">pvleuven</a> helped fixing a warning in ICC.</li>
<li><a href="https://github.com/crea7or">Pavel</a> helped fixing some warnings in MSVC.</li>
<li><a href="https://github.com/jseward">Jamie Seward</a> avoided unnecessary string copies in <code>find()</code> and <code>count()</code>.</li>
<li><a href="https://github.com/Itja">Mitja</a> fixed some typos.</li>
<li><a href="https://github.com/jowr">Jorrit Wronski</a> updated the Hunter package links.</li>
<li><a href="https://github.com/TinyTinni">Matthias Möller</a> added a <code>.natvis</code> for the MSVC debug view.</li>
<li><a href="https://github.com/bogemic">bogemic</a> fixed some C++17 deprecation warnings.</li>
<li><a href="https://github.com/erengy">Eren Okka</a> fixed some MSVC warnings.</li>
<li><a href="https://github.com/abolz">abolz</a> integrated the Grisu2 algorithm for proper floating-point formatting, allowing more roundtrip checks to succeed.</li>
<li><a href="https://github.com/Pipeliner">Vadim Evard</a> fixed a Markdown issue in the README.</li>
<li><a href="https://github.com/zerodefect">zerodefect</a> fixed a compiler warning.</li>
<li><a href="https://github.com/kaidokert">Kert</a> allowed to template the string type in the serialization and added the possibility to override the exceptional behavior.</li>
<li><a href="https://github.com/mark-99">mark-99</a> helped fixing an ICC error.</li>
<li><a href="https://github.com/patrikhuber">Patrik Huber</a> fixed links in the README file.</li>
<li><a href="https://github.com/johnfb">johnfb</a> found a bug in the implementation of CBOR's indefinite length strings.</li>
<li><a href="https://github.com/pfultz2">Paul Fultz II</a> added a note on the cget package manager.</li>
<li><a href="https://github.com/wla80">Wilson Lin</a> made the integration section of the README more concise.</li>
<li><a href="https://github.com/ralfbielig">RalfBielig</a> detected and fixed a memory leak in the parser callback.</li>
<li><a href="https://github.com/agrianius">agrianius</a> allowed to dump JSON to an alternative string type.</li>
<li><a href="https://github.com/ktonon">Kevin Tonon</a> overworked the C++11 compiler checks in CMake.</li>
<li><a href="https://github.com/ax3l">Axel Huebl</a> simplified a CMake check and added support for the <a href="https://spack.io" rel="nofollow">Spack package manager</a>.</li>
<li><a href="https://github.com/coryan">Carlos O'Ryan</a> fixed a typo.</li>
<li><a href="https://github.com/jammehcow">James Upjohn</a> fixed a version number in the compilers section.</li>
<li><a href="https://github.com/chuckatkins">Chuck Atkins</a> adjusted the CMake files to the CMake packaging guidelines and provided documentation for the CMake integration.</li>
<li><a href="https://github.com/dns13">Jan Schöppach</a> fixed a typo.</li>
<li><a href="https://github.com/martin-mfg">martin-mfg</a> fixed a typo.</li>
<li><a href="https://github.com/TinyTinni">Matthias Möller</a> removed the dependency from <code>std::stringstream</code>.</li>
<li><a href="https://github.com/agrianius">agrianius</a> added code to use alternative string implementations.</li>
<li><a href="https://github.com/Daniel599">Daniel599</a> allowed to use more algorithms with the <code>items()</code> function.</li>
<li><a href="https://github.com/jrakow">Julius Rakow</a> fixed the Meson include directory and fixed the links to <a href="/nlohmann/json/blob/develop/cppreference.com">cppreference.com</a>.</li>
<li><a href="https://github.com/sonulohani">Sonu Lohani</a> fixed the compilation with MSVC 2015 in debug mode.</li>
<li><a href="https://github.com/grembo">grembo</a> fixed the test suite and re-enabled several test cases.</li>
<li><a href="https://github.com/simnalamburt">Hyeon Kim</a> introduced the macro <code>JSON_INTERNAL_CATCH</code> to control the exception handling inside the library.</li>
<li><a href="https://github.com/thyu">thyu</a> fixed a compiler warning.</li>
<li><a href="https://github.com/LEgregius">David Guthrie</a> fixed a subtle compilation error with Clang 3.4.2.</li>
<li><a href="https://github.com/dennisfischer">Dennis Fischer</a> allowed to call <code>find_package</code> without installing the library.</li>
<li><a href="https://github.com/simnalamburt">Hyeon Kim</a> fixed an issue with a double macro definition.</li>
<li><a href="https://github.com/rivertam">Ben Berman</a> made some error messages more understandable.</li>
<li><a href="https://github.com/zakalibit">zakalibit</a> fixed a compilation problem with the Intel C++ compiler.</li>
<li><a href="https://github.com/mandreyel">mandreyel</a> fixed a compilation problem.</li>
<li><a href="https://github.com/koponomarenko">Kostiantyn Ponomarenko</a> added version and license information to the Meson build file.</li>
<li><a href="https://github.com/henryiii">Henry Schreiner</a> added support for GCC 4.8.</li>
<li><a href="https://github.com/knilch0r">knilch</a> made sure the test suite does not stall when run in the wrong directory.</li>
<li><a href="https://github.com/antonioborondo">Antonio Borondo</a> fixed an MSVC 2017 warning.</li>
<li><a href="https://github.com/dgendreau">Dan Gendreau</a> implemented the <code>NLOHMANN_JSON_SERIALIZE_ENUM</code> macro to quickly define a enum/JSON mapping.</li>
<li><a href="https://github.com/efp">efp</a> added line and column information to parse errors.</li>
<li><a href="https://github.com/julian-becker">julian-becker</a> added BSON support.</li>
</ul>
<p>Thanks a lot for helping out! Please <a href="mailto:mail@nlohmann.me">let me know</a> if I forgot someone.</p>
<h2><a id="user-content-used-third-party-tools" class="anchor" aria-hidden="true" href="#used-third-party-tools"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Used third-party tools</h2>
<p>The library itself consists of a single header file licensed under the MIT license. However, it is built, tested, documented, and whatnot using a lot of third-party tools and services. Thanks a lot!</p>
<ul>
<li><a href="https://github.com/edlund/amalgamate"><strong>amalgamate.py - Amalgamate C source and header files</strong></a> to create a single header file</li>
<li><a href="http://lcamtuf.coredump.cx/afl/" rel="nofollow"><strong>American fuzzy lop</strong></a> for fuzz testing</li>
<li><a href="https://www.appveyor.com" rel="nofollow"><strong>AppVeyor</strong></a> for <a href="https://ci.appveyor.com/project/nlohmann/json" rel="nofollow">continuous integration</a> on Windows</li>
<li><a href="http://astyle.sourceforge.net" rel="nofollow"><strong>Artistic Style</strong></a> for automatic source code identation</li>
<li><a href="https://github.com/philsquared/Catch"><strong>Catch</strong></a> for the unit tests</li>
<li><a href="http://clang.llvm.org" rel="nofollow"><strong>Clang</strong></a> for compilation with code sanitizers</li>
<li><a href="https://cmake.org" rel="nofollow"><strong>CMake</strong></a> for build automation</li>
<li><a href="https://www.codacy.com" rel="nofollow"><strong>Codacity</strong></a> for further <a href="https://www.codacy.com/app/nlohmann/json" rel="nofollow">code analysis</a></li>
<li><a href="https://coveralls.io" rel="nofollow"><strong>Coveralls</strong></a> to measure <a href="https://coveralls.io/github/nlohmann/json" rel="nofollow">code coverage</a></li>
<li><a href="https://scan.coverity.com" rel="nofollow"><strong>Coverity Scan</strong></a> for <a href="https://scan.coverity.com/projects/nlohmann-json" rel="nofollow">static analysis</a></li>
<li><a href="http://cppcheck.sourceforge.net" rel="nofollow"><strong>cppcheck</strong></a> for static analysis</li>
<li><a href="http://www.stack.nl/~dimitri/doxygen/" rel="nofollow"><strong>Doxygen</strong></a> to generate <a href="https://nlohmann.github.io/json/" rel="nofollow">documentation</a></li>
<li><a href="https://github.com/rstacruz/git-update-ghpages"><strong>git-update-ghpages</strong></a> to upload the documentation to gh-pages</li>
<li><a href="https://github.com/skywinder/github-changelog-generator"><strong>GitHub Changelog Generator</strong></a> to generate the <a href="https://github.com/nlohmann/json/blob/develop/ChangeLog.md">ChangeLog</a></li>
<li><a href="https://github.com/google/benchmark"><strong>Google Benchmark</strong></a> to implement the benchmarks</li>
<li><a href="http://llvm.org/docs/LibFuzzer.html" rel="nofollow"><strong>libFuzzer</strong></a> to implement fuzz testing for OSS-Fuzz</li>
<li><a href="https://github.com/google/oss-fuzz"><strong>OSS-Fuzz</strong></a> for continuous fuzz testing of the library (<a href="https://github.com/google/oss-fuzz/tree/master/projects/json">project repository</a>)</li>
<li><a href="https://probot.github.io" rel="nofollow"><strong>Probot</strong></a> for automating maintainer tasks such as closing stale issues, requesting missing information, or detecting toxic comments.</li>
<li><a href="https://github.com/nlohmann/json/blob/develop/doc/scripts/send_to_wandbox.py"><strong>send_to_wandbox</strong></a> to send code examples to <a href="http://melpon.org/wandbox" rel="nofollow">Wandbox</a></li>
<li><a href="https://travis-ci.org" rel="nofollow"><strong>Travis</strong></a> for <a href="https://travis-ci.org/nlohmann/json" rel="nofollow">continuous integration</a> on Linux and macOS</li>
<li><a href="http://valgrind.org" rel="nofollow"><strong>Valgrind</strong></a> to check for correct memory management</li>
<li><a href="http://melpon.org/wandbox" rel="nofollow"><strong>Wandbox</strong></a> for <a href="https://wandbox.org/permlink/TarF5pPn9NtHQjhf" rel="nofollow">online examples</a></li>
</ul>
<h2><a id="user-content-projects-using-json-for-modern-c" class="anchor" aria-hidden="true" href="#projects-using-json-for-modern-c"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Projects using JSON for Modern C++</h2>
<p>The library is currently used in Apple macOS Sierra and iOS 10. I am not sure what they are using the library for, but I am happy that it runs on so many devices.</p>
<h2><a id="user-content-notes" class="anchor" aria-hidden="true" href="#notes"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Notes</h2>
<ul>
<li>The code contains numerous debug <strong>assertions</strong> which can be switched off by defining the preprocessor macro <code>NDEBUG</code>, see the <a href="https://en.cppreference.com/w/cpp/error/assert" rel="nofollow">documentation of <code>assert</code></a>. In particular, note <a href="https://nlohmann.github.io/json/classnlohmann_1_1basic__json_a2e26bd0b0168abb61f67ad5bcd5b9fa1.html#a2e26bd0b0168abb61f67ad5bcd5b9fa1" rel="nofollow"><code>operator[]</code></a> implements <strong>unchecked access</strong> for const objects: If the given key is not present, the behavior is undefined (think of a dereferenced null pointer) and yields an <a href="https://github.com/nlohmann/json/issues/289">assertion failure</a> if assertions are switched on. If you are not sure whether an element in an object exists, use checked access with the <a href="https://nlohmann.github.io/json/classnlohmann_1_1basic__json_a674de1ee73e6bf4843fc5dc1351fb726.html#a674de1ee73e6bf4843fc5dc1351fb726" rel="nofollow"><code>at()</code> function</a>.</li>
<li>As the exact type of a number is not defined in the <a href="http://rfc7159.net/rfc7159" rel="nofollow">JSON specification</a>, this library tries to choose the best fitting C++ number type automatically. As a result, the type <code>double</code> may be used to store numbers which may yield <a href="https://github.com/nlohmann/json/issues/181"><strong>floating-point exceptions</strong></a> in certain rare situations if floating-point exceptions have been unmasked in the calling code. These exceptions are not caused by the library and need to be fixed in the calling code, such as by re-masking the exceptions prior to calling library functions.</li>
<li>The library supports <strong>Unicode input</strong> as follows:
<ul>
<li>Only <strong>UTF-8</strong> encoded input is supported which is the default encoding for JSON according to <a href="http://rfc7159.net/rfc7159#rfc.section.8.1" rel="nofollow">RFC 7159</a>.</li>
<li>Other encodings such as Latin-1, UTF-16, or UTF-32 are not supported and will yield parse or serialization errors.</li>
<li><a href="http://www.unicode.org/faq/private_use.html#nonchar1" rel="nofollow">Unicode noncharacters</a> will not be replaced by the library.</li>
<li>Invalid surrogates (e.g., incomplete pairs such as <code>\uDEAD</code>) will yield parse errors.</li>
<li>The strings stored in the library are UTF-8 encoded. When using the default string type (<code>std::string</code>), note that its length/size functions return the number of stored bytes rather than the number of characters or glyphs.</li>
</ul>
</li>
<li>The code can be compiled without C++ <strong>runtime type identification</strong> features; that is, you can use the <code>-fno-rtti</code> compiler flag.</li>
<li><strong>Exceptions</strong> are used widely within the library. They can, however, be switched off with either using the compiler flag <code>-fno-exceptions</code> or by defining the symbol <code>JSON_NOEXCEPTION</code>. In this case, exceptions are replaced by an <code>abort()</code> call.</li>
<li>By default, the library does not preserve the <strong>insertion order of object elements</strong>. This is standards-compliant, as the <a href="https://tools.ietf.org/html/rfc7159.html" rel="nofollow">JSON standard</a> defines objects as "an unordered collection of zero or more name/value pairs". If you do want to preserve the insertion order, you can specialize the object type with containers like <a href="https://github.com/Tessil/ordered-map"><code>tsl::ordered_map</code></a> (<a href="https://github.com/nlohmann/json/issues/546#issuecomment-304447518">integration</a>) or <a href="https://github.com/nlohmann/fifo_map"><code>nlohmann::fifo_map</code></a> (<a href="https://github.com/nlohmann/json/issues/485#issuecomment-333652309">integration</a>).</li>
</ul>
<h2><a id="user-content-execute-unit-tests" class="anchor" aria-hidden="true" href="#execute-unit-tests"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9h1v1H4c-1.5 0-3-1.69-3-3.5S2.55 3 4 3h4c1.45 0 3 1.69 3 3.5 0 1.41-.91 2.72-2 3.25V8.59c.58-.45 1-1.27 1-2.09C10 5.22 8.98 4 8 4H4c-.98 0-2 1.22-2 2.5S3 9 4 9zm9-3h-1v1h1c1 0 2 1.22 2 2.5S13.98 12 13 12H9c-.98 0-2-1.22-2-2.5 0-.83.42-1.64 1-2.09V6.25c-1.09.53-2 1.84-2 3.25C6 11.31 7.55 13 9 13h4c1.45 0 3-1.69 3-3.5S14.5 6 13 6z"></path></svg></a>Execute unit tests</h2>
<p>To compile and run the tests, you need to execute</p>
<div class="highlight highlight-source-shell"><pre>$ mkdir build
$ <span class="pl-c1">cd</span> build
$ cmake ..
$ cmake --build <span class="pl-c1">.</span>
$ ctest --output-on-failure</pre></div>
<p>For more information, have a look at the file <a href="https://github.com/nlohmann/json/blob/master/.travis.yml">.travis.yml</a>.</p>
</article>
  </div>

    </div>

  

  <details class="details-reset details-overlay details-overlay-dark">
    <summary data-hotkey="l" aria-label="Jump to line"></summary>
    <details-dialog class="Box Box--overlay d-flex flex-column anim-fade-in fast linejump" aria-label="Jump to line">
      <!-- '"` --><!-- </textarea></xmp> --></option></form><form class="js-jump-to-line-form Box-body d-flex" action="" accept-charset="UTF-8" method="get"><input name="utf8" type="hidden" value="&#x2713;" />
        <input class="form-control flex-auto mr-3 linejump-input js-jump-to-line-field" type="text" placeholder="Jump to line&hellip;" aria-label="Jump to line" autofocus>
        <button type="submit" class="btn" data-close-dialog>Go</button>
</form>    </details-dialog>
  </details>


  </div>
  <div class="modal-backdrop js-touch-events"></div>
</div>

    </div>
  </div>

  </div>

        
<div class="footer container-lg px-3" role="contentinfo">
  <div class="position-relative d-flex flex-justify-between pt-6 pb-2 mt-6 f6 text-gray border-top border-gray-light ">
    <ul class="list-style-none d-flex flex-wrap ">
      <li class="mr-3">&copy; 2018 <span title="0.10915s from unicorn-7bd959dcd-dczjg">GitHub</span>, Inc.</li>
        <li class="mr-3"><a data-ga-click="Footer, go to terms, text:terms" href="https://github.com/site/terms">Terms</a></li>
        <li class="mr-3"><a data-ga-click="Footer, go to privacy, text:privacy" href="https://github.com/site/privacy">Privacy</a></li>
        <li class="mr-3"><a href="https://help.github.com/articles/github-security/" data-ga-click="Footer, go to security, text:security">Security</a></li>
        <li class="mr-3"><a href="https://status.github.com/" data-ga-click="Footer, go to status, text:status">Status</a></li>
        <li><a data-ga-click="Footer, go to help, text:help" href="https://help.github.com">Help</a></li>
    </ul>

    <a aria-label="Homepage" title="GitHub" class="footer-octicon mr-lg-4" href="https://github.com">
      <svg height="24" class="octicon octicon-mark-github" viewBox="0 0 16 16" version="1.1" width="24" aria-hidden="true"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"/></svg>
</a>
   <ul class="list-style-none d-flex flex-wrap ">
        <li class="mr-3"><a data-ga-click="Footer, go to contact, text:contact" href="https://github.com/contact">Contact GitHub</a></li>
        <li class="mr-3"><a href="https://github.com/pricing" data-ga-click="Footer, go to Pricing, text:Pricing">Pricing</a></li>
      <li class="mr-3"><a href="https://developer.github.com" data-ga-click="Footer, go to api, text:api">API</a></li>
      <li class="mr-3"><a href="https://training.github.com" data-ga-click="Footer, go to training, text:training">Training</a></li>
        <li class="mr-3"><a href="https://blog.github.com" data-ga-click="Footer, go to blog, text:blog">Blog</a></li>
        <li><a data-ga-click="Footer, go to about, text:about" href="https://github.com/about">About</a></li>

    </ul>
  </div>
  <div class="d-flex flex-justify-center pb-6">
    <span class="f6 text-gray-light"></span>
  </div>
</div>



  <div id="ajax-error-message" class="ajax-error-message flash flash-error">
    <svg class="octicon octicon-alert" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M8.893 1.5c-.183-.31-.52-.5-.887-.5s-.703.19-.886.5L.138 13.499a.98.98 0 0 0 0 1.001c.193.31.53.501.886.501h13.964c.367 0 .704-.19.877-.5a1.03 1.03 0 0 0 .01-1.002L8.893 1.5zm.133 11.497H6.987v-2.003h2.039v2.003zm0-3.004H6.987V5.987h2.039v4.006z"/></svg>
    <button type="button" class="flash-close js-ajax-error-dismiss" aria-label="Dismiss error">
      <svg class="octicon octicon-x" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48L7.48 8z"/></svg>
    </button>
    You can’t perform that action at this time.
  </div>


    <script crossorigin="anonymous" integrity="sha512-WnyO4VoIUwWWQOmFLjYf4UGg/c1z9VlaLN8IMuiI3uMhhl6rejyThRdLPDyePeUPW6N+38OoBMs6AkqcvWALtA==" type="application/javascript" src="https://assets-cdn.github.com/assets/compat-b66b5d97b4442a01f057c74b091c4368.js"></script>
    <script crossorigin="anonymous" integrity="sha512-/iZWaUnKp78VsMWZGo2aU/WzWHTwFbZrhAGxpWOtkLb/dRMR14ef1rQsgmEK5Lfvk9AtEUidXFgLPfAnczjI2g==" type="application/javascript" src="https://assets-cdn.github.com/assets/frameworks-f15af08e5888bd0158d6972a2ff991e9.js"></script>
    
    <script crossorigin="anonymous" async="async" integrity="sha512-5h5Hx6SK1gdmCM3D+3XTE7ROlh5ZAQbZkQvRmWGlCoLrcWH6C9Cv4fyt0JH3ZHJKDPk6orSchC3TArxh0uuCZg==" type="application/javascript" src="https://assets-cdn.github.com/assets/github-a2c56bc2a974067bfd56642aa2108f23.js"></script>
    
    
    
  <div class="js-stale-session-flash stale-session-flash flash flash-warn flash-banner d-none">
    <svg class="octicon octicon-alert" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M8.893 1.5c-.183-.31-.52-.5-.887-.5s-.703.19-.886.5L.138 13.499a.98.98 0 0 0 0 1.001c.193.31.53.501.886.501h13.964c.367 0 .704-.19.877-.5a1.03 1.03 0 0 0 .01-1.002L8.893 1.5zm.133 11.497H6.987v-2.003h2.039v2.003zm0-3.004H6.987V5.987h2.039v4.006z"/></svg>
    <span class="signed-in-tab-flash">You signed in with another tab or window. <a href="">Reload</a> to refresh your session.</span>
    <span class="signed-out-tab-flash">You signed out in another tab or window. <a href="">Reload</a> to refresh your session.</span>
  </div>
  <div class="facebox" id="facebox" style="display:none;">
  <div class="facebox-popup">
    <div class="facebox-content" role="dialog" aria-labelledby="facebox-header" aria-describedby="facebox-description">
    </div>
    <button type="button" class="facebox-close js-facebox-close" aria-label="Close modal">
      <svg class="octicon octicon-x" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48L7.48 8z"/></svg>
    </button>
  </div>
</div>

  <template id="site-details-dialog">
  <details class="details-reset details-overlay details-overlay-dark lh-default text-gray-dark" open>
    <summary aria-haspopup="dialog" aria-label="Close dialog"></summary>
    <details-dialog class="Box Box--overlay d-flex flex-column anim-fade-in fast">
      <button class="Box-btn-octicon m-0 btn-octicon position-absolute right-0 top-0" type="button" aria-label="Close dialog" data-close-dialog>
        <svg class="octicon octicon-x" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48L7.48 8z"/></svg>
      </button>
      <div class="octocat-spinner my-6 js-details-dialog-spinner"></div>
    </details-dialog>
  </details>
</template>

  <div class="Popover js-hovercard-content position-absolute" style="display: none; outline: none;" tabindex="0">
  <div class="Popover-message Popover-message--bottom-left Popover-message--large Box box-shadow-large" style="width:360px;">
  </div>
</div>

<div id="hovercard-aria-description" class="sr-only">
  Press h to open a hovercard with more details.
</div>


  </body>
</html>

