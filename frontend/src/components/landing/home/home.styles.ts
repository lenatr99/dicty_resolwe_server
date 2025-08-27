import styled from 'styled-components';
import { Button } from '@mui/material';
import dictyImage from 'images/dicty_image.png';
import { breakpoints } from 'components/app/globalStyle';

export const HomeContainer = styled.div`
    background: url(${dictyImage}), linear-gradient(to right, #49688d, #7ca2cf);
    background-repeat: no-repeat;
    background-position: left bottom;
    padding-bottom: 70px;
    width: 100%;
    position: relative;
    color: #fff;
`;

export const HomeWrapper = styled.div`
    width: 100%;
    display: flex;
    flex-flow: column nowrap;

    @media (min-width: ${breakpoints.small}px) {
        flex-flow: row wrap;
        justify-content: space-between;
    }
`;

export const HomeDescription = styled.div`
    width: 100%;

    & p {
        font-size: 1.75rem;
        line-height: 1.2em;
        margin: 1rem 0 3rem;
    }

    > .title {
        display: flex;
        align-items: baseline;
        justify-content: flex-start;
        flex-wrap: wrap;
        margin-top: -1.7rem;
        margin-bottom: 2rem;

        > * {
            margin: 0 0.5rem 0.5rem 0;
        }

        .logo {
            position: relative;
            top: 1.7rem;
            flex-wrap: nowrap;
        }

        .header {
            margin-top: 1rem;
            margin-bottom: 0;
        }

        .version {
            white-space: nowrap;
            line-height: 2rem;
        }
    }

    @media (min-width: ${breakpoints.small}px) {
        width: 59%;
    }
`;

export const HomeDemo = styled.div`
    width: 100%;
    display: flex;
    flex-direction: column;
    justify-content: flex-start;

    & > img {
        max-width: 100%;
        height: auto;
        padding-bottom: 1rem;
    }

    @media (min-width: ${breakpoints.small}px) {
        width: 40%;
    }
`;

export const HomeRunDictyButton = styled(Button)`
    color: #fff;
    border-color: #fff;
    background-color: #49688d;
    background-image: none;
    text-shadow: none;
    font-size: 30px;
    text-transform: none;

    & a {
        text-decoration: none;
        color: inherit;
    }
`;
