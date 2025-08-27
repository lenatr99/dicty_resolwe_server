import styled from 'styled-components';
import { Link } from '@mui/material';

export const FooterContainer = styled.div`
    background-color: #6f6f6f;
    padding-top: 5px;
    width: 100%;
    color: #fff;

    & a {
        color: #afafaf;
    }

    & h2,
    h3,
    h4,
    h5,
    h6 {
        font-family: 'FS Joey Web Bold', Helvetica, Arial, Verdana, sans-serif;
    }

    & p {
        font-size: 1.08em;
    }
`;

export const FooterWrapper = styled.div`
    width: 100%;
    padding-top: 35px;
    padding-bottom: 30px;
    font-size: 1.143em;
    display: flex;
    flex-flow: row wrap;
    justify-content: space-between;
    align-items: center;
`;

export const SocialNetworkLink = styled(Link)`
    margin-right: 10px;

    & svg {
        font-size: 48px;
    }
    & svg:hover {
        color: #fff;
    }
`;

export const SocialNetworkLinkWrapper = styled.span`
    margin-right: 10px;
`;

export const FooterTitle = styled.h2`
    margin: 0;
`;
