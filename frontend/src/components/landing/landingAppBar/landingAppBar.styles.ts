import styled from 'styled-components';
import { Button, AppBar } from '@mui/material';

export const LandingAppBarWrapper = styled(AppBar)`
    background-color: #fff;
`;

export const ShortDescription = styled.div`
    margin-left: auto;
    color: #49688d;
    display: inline-flex;
    white-space: pre;
`;

export const NavButton = styled(Button)`
    text-transform: none;
    justify-content: flex-start;
`;
