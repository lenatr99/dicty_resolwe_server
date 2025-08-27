import styled from 'styled-components';
import { breakpoints } from 'components/app/globalStyle';

export const CreditsLogos = styled.div`
    width: 100%;
    display: flex;
    flex-flow: column nowrap;
    margin-top: 50px;

    & img {
        margin-bottom: 40px;
    }

    @media (min-width: ${breakpoints.small}px) {
        flex-flow: row wrap;
        justify-content: center;

        & img {
            margin-right: 35px;
            margin-left: 35px;
        }
    }
`;
