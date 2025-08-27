import styled from 'styled-components';
import { breakpoints } from 'components/app/globalStyle';

export const ReferencesContainer = styled.div`
    width: 100%;
    display: flex;
    flex-flow: column nowrap;

    @media (min-width: 768px) {
        flex-flow: row wrap;
        justify-content: space-between;
    }
`;

export const ReferenceContainer = styled.div`
    width: 100%;
    text-align: center;
    display: flex;
    flex-direction: column;
    align-content: center;
    margin-bottom: 35px;

    @media (min-width: ${breakpoints.small}px) {
        width: 49%;
    }
`;

export const ReferenceAuthors = styled.p`
    font-size: 1.1em;
`;

export const ReferenceSource = styled.p`
    color: #333;
`;
