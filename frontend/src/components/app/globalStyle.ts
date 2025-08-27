import { createGlobalStyle } from 'styled-components';

export const GlobalStyle = createGlobalStyle`
  body {
    margin: 0;
    min-width: 360px;
    padding-top: 0px;
    font-family: "FS Joey Web Regular",Helvetica,Arial,Verdana,sans-serif;
    font-size: 1rem;
    font-style: normal;
    font-weight: normal;
  }

  h1 {
    font-family: "FS Joey Web Bold",Helvetica,Arial,Verdana,sans-serif;
    font-size: 4rem;
    margin-top: 10px;
    margin-bottom:20px;
  }

  h2, h3, h4, h5, h6 {
    font-family: "FS Joey Web Regular",Helvetica,Arial,Verdana,sans-serif;
    margin-top: 20px;
    margin-bottom: 10px;
  }

  p {
    font-size: 1.3rem;
    line-height: 1.2em;
    margin: 0 0 10px;
  }

  a {
    color: #428bca;
    text-decoration: none;
  }
`;

export const breakpoints = {
    small: 768,
    mid: 992,
    large: 1200,
};

export const appBarHeight = 56;
