import layoutsReducer, {
    generateModuleLayout,
    layoutsChanged,
    layoutsInitialState,
    LayoutsState,
} from './layouts';
import { LayoutBreakpoint, ModulesKeys } from 'components/genexpress/common/constants';
import { BreakpointsCols } from 'redux/models/internal';

const widthInCols = 3;
const minWidthInCols = 2;
const heightInCols = 3;
const minHeightInCols = 2;
const breakpointCols: BreakpointsCols = { large: 12, mid: 10, small: 6 };
const moduleKey = ModulesKeys.expressionTimeCourses;

describe('layouts store', () => {
    let initialState: LayoutsState;

    beforeEach(() => {
        initialState = {
            layouts: [],
        };
    });

    it('should change state with layoutsChanged action', () => {
        const newState = layoutsReducer(initialState, layoutsChanged(layoutsInitialState));

        expect(newState).toEqual(layoutsInitialState);
    });

    it('should put the module in same row if space is available', () => {
        const availableSpace = 5;
        const lastRowAvailableSpaceInCols = {
            rowIndex: 1,
            availableSpaceInCols: availableSpace,
        };

        const moduleLayout = generateModuleLayout(
            moduleKey,
            LayoutBreakpoint.large,
            lastRowAvailableSpaceInCols,
            widthInCols,
            minWidthInCols,
            heightInCols,
            minHeightInCols,
            breakpointCols,
        );

        expect(moduleLayout).toEqual({
            i: moduleKey,
            x: 7,
            y: heightInCols,
            w: widthInCols,
            h: 1 * heightInCols, // 1 === rowIndex
            minW: minWidthInCols,
            minH: minHeightInCols,
        });

        expect(lastRowAvailableSpaceInCols).toEqual({
            rowIndex: 1,
            availableSpaceInCols: availableSpace - widthInCols,
        });
    });

    it('should put the module in the next row if space is not available', () => {
        const lastRowAvailableSpaceInCols = {
            rowIndex: 1,
            availableSpaceInCols: 2,
        };

        const moduleLayout = generateModuleLayout(
            moduleKey,
            LayoutBreakpoint.large,
            lastRowAvailableSpaceInCols,
            widthInCols,
            minWidthInCols,
            heightInCols,
            minHeightInCols,
            breakpointCols,
        );

        expect(moduleLayout).toEqual({
            i: moduleKey,
            x: 0,
            y: 2 * heightInCols, // 2 === rowIndex
            w: widthInCols,
            h: heightInCols,
            minW: minWidthInCols,
            minH: minHeightInCols,
        });

        expect(lastRowAvailableSpaceInCols).toEqual({
            rowIndex: 2,
            availableSpaceInCols: breakpointCols.large - widthInCols,
        });
    });
});
